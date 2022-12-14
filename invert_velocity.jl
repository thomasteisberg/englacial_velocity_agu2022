using HDF5
using Plots
using DifferentialEquations
using Dierckx
using Statistics
using ColorSchemes

#plotlyjs() # Nice for interactive plots to explore
pythonplot() # Used for figure generation

#data = h5open("model_output_sim_n3.5_fric100_Btemp-10.mat", "r")["output_data"] # works fine, essentially plug flow
#data = h5open("model_output_sim_n2.5_fric200_Btemp-20.mat", "r")["output_data"] # unstable, but not clear why

# Included example data files:
data = h5open("issm/sim_data_n3.5_fric100_Btemp-10.mat", "r")["output_data"] # Essentially plug flow
#data = h5open("issm/sim_data_n2.5_fric200_Btemp-20.mat", "r")["output_data"] # Shows more variation with depth


# Estimate from data:
# -- surface velocity --
# Vxs, Vys -- surface velocity
# dVxs/dx, dVys/dy -- spatial derivatives of surfave velocity
# -- layer velocity --
# (l1,l2)_Vz -- vertical velocity of layer above and below
# L -- layer thickness
# dL/dt -- change in layer thickness over time
# dL/dx -- layer slope

# NOTE: The extracted slice from the model output flows along the Y axis.
# In the ODE, we use X as the along-flow spatial dimension.
# This just means that we swap Y for X in the code here.


# Extract ys vector and convert it to a Julia range
ys_data = data["layers"]["l0"]["ys"][1,:]
ys_spacing = mean(diff(ys_data))
@assert maximum(diff(ys_data)) < (ys_spacing + 1e-3)
@assert minimum(diff(ys_data)) > (ys_spacing - 1e-3)
ys = minimum(ys_data):ys_spacing:maximum(ys_data)
@assert length(ys_data) == length(ys)

# Surface velocity
Vxs = data["layers"]["l0"]["vx"][1,:]
Vys = data["layers"]["l0"]["vy"][1,:]
#   interpolators
Vys_interp = Spline1D(ys, Vys, s=1e-1)
Vxs_interp = Spline1D(ys, Vxs)

# Surface velocity spatial derivatives
#   dVys/dy comes directly from the interpolator
dVys_dy_grad(y) = derivative(Vys_interp, y)
dVys_dy_interp(ys) = [dVys_dy_grad(y)[1] for y in ys]
#   dVxs/dx has to be estimated by finite differences
x_slice_pos_x = data["slices"]["offset_pos"]["x_slice"][1,1]
x_slice_neg_x = data["slices"]["offset_neg"]["x_slice"][1,1]
slices_delta_x = x_slice_pos_x - x_slice_neg_x
delta_Vxs = data["layers"]["l0"]["vx_offset_pos"][1,:] - data["layers"]["l0"]["vx_offset_neg"][1,:]
dVxs_dx_fd = delta_Vxs / slices_delta_x
mask = .!isnan.(dVxs_dx_fd)
dVxs_dx_interp = Spline1D(ys[mask], dVxs_dx_fd[mask], s=1e-1)

# Surface velocity plot
# begin
#     #plot(ys, Vys, label="Vys")
#     plot(ys, dVys_dy_interp(ys), label="dVys_dy")
#     plot!(ys, dVxs_dx_interp(ys), label="dVxs_dx")
#     xlims!(0, 20e3)
#     #ylims!(-0.5e-3, 1e-3)
# end

# Solve for horizontal velocity between two layers
function solve_between_layers(selected_layers)
    l1 = data["layers"][selected_layers[1]]
    l2 = data["layers"][selected_layers[2]]
    
    l1_z_interp = Spline1D(l1["ys"][1,:], l1["zs"][1,:])
    l2_z_interp = Spline1D(l2["ys"][1,:], l2["zs"][1,:])

    # Create mapping from w (path length along midpoint between layers) to y
    #midpoint_z_interp(y) = (l2_z_interp(y) + l1_z_interp(y)) / 2
    # Derivative of z value along the path w/r to the horizontal coordinate
    path_dz_dx(y) = (derivative(l2_z_interp, y) + derivative(l1_z_interp, y)) / 2


    l1_vz = Spline1D(l1["ys"][1,:], l1["vz"][1,:])
    l2_vz = Spline1D(l2["ys"][1,:], l2["vz"][1,:])
    L_interp = Spline1D(l1["ys"][1,:], l2_z_interp(l1["ys"][1,:]) .- l1_z_interp(l1["ys"][1,:]))
    dL_dx_grad(y) = derivative(L_interp, y)
    dL_dx_interp(ys) = [dL_dx_grad(y)[1] for y in ys]
    dL_dt_interp = Spline1D(ys, (l2_vz(ys) .- l1_vz(ys)))
    # Note: above should be divided by time period between observations,
    # but, in this case, that is 1 year

    begin # Generate some debugging plots
        ws = 0:100:25e3
        p1 = plot(ws, dVxs_dx_interp(ws), label="dVxs_dx")
        plot!(ws, dVys_dy_interp(ws), label="dVys_dy")

        plot!(ws, l1_vz(ws), label="l1_vz")
        plot!(ws, l2_vz(ws), label="l2_vz")
        
        plot!(ws, dL_dt_interp(ws), label="dL_dt")
        plot!(ws, dL_dx_interp(ws), label="dL_dx")

        p2 = plot(ws, l1_z_interp(ws), label="l1 z")
        plot!(ws, l2_z_interp(ws), label="l2 z")
        #plot!(ws, path_dz_dx(ws), label="path_dz_dx")
        plot!(ws, L_interp(ws), label="L")

        p3 = plot(ws, Vys_interp(ws), label="Vys")

        stab_crit_b_thick_part = (dL_dx_interp(ws) ./ L_interp(ws))
        stab_crit_b_div_part = ((dVxs_dx_interp(ws) + dVys_dy_interp(ws)) ./ Vys_interp(ws))
        stab_crit_b = stab_crit_b_thick_part + stab_crit_b_div_part
        
        p4 = plot(ws, stab_crit_b, label="b")
        plot!(ws, stab_crit_b_thick_part, label="dL/dx * 1/L")
        plot!(ws, stab_crit_b_div_part, label="(dVxs_dx+dVys_dy)/Vys")
        plot!(ws, dVys_dy_interp(ws), label="dVys_dy")
        plot!(ws, dVxs_dx_interp(ws), label="dVxs_dx")

        a = (1 ./ Vys_interp(ws)) .* (1 ./ L_interp(ws)) .* dL_dt_interp(ws)
        p5 = plot(ws, a, label="a")
        plot!(ws, (1 ./ Vys_interp(ws)), label="1/Vys")
        plot!(ws, (1 ./ L_interp(ws)), label="1/L")
        plot!(ws, dL_dt_interp(ws), label="dL/dt")
    end

    # Define ODE problem
    function ds_dw(s, p, w) # u, p t
        dVxs_dx = dVys_dy_interp(w)[1]
        dVys_dy = dVxs_dx_interp(w)[1]

        Vxs = Vys_interp(w)
        L = L_interp(w)
        dL_dt = dL_dt_interp(w)
        dL_dx = dL_dx_interp(w)[1]

        dz_dx = path_dz_dx(w)[1]

        # This is the scalar field value that we're evaluating
        scalar_field_value = (-1 / Vxs) * ( (1/L) * (dL_dt + s * Vxs * dL_dx) + s * (dVxs_dx + dVys_dy))
        # Solving the ODE on a curved path, so we need to account for path length
        # See https://en.wikipedia.org/wiki/Line_integral#Line_integral_of_a_scalar_field
        path_derivative = sqrt(1 + dz_dx)

        return scalar_field_value * path_derivative
    end

    s0 = 0.5
    wspan = (10, 20e3)
    prob = ODEProblem(ds_dw, s0, wspan)
    sol = solve(prob)

    results = Dict(
        "selected_layers" => selected_layers,
        "l1" => l1,
        "l2" => l2,
        "l1_z" => l1_z_interp,
        "l2_z" => l2_z_interp,
        "solution" => sol,
        "plots" => [p1, p2, p3, p4, p5]
    )
    return results
end

begin # Compute layer velocity results for each layer pair
    layer_isochrone_year(layer_string) = parse(Int, layer_string[2:end])

    layer_keys = sort(keys(data["layers"]), by=layer_isochrone_year)
    layer_keys = layer_keys[2:end]
    
    results = Vector{Dict}()
    for k in 1:length(layer_keys)-1
        selected_layers = [layer_keys[k+1], layer_keys[k]]
        push!(results, solve_between_layers(selected_layers))
        println("Layers: $(selected_layers), \tSolution result: $(results[k]["solution"].retcode)")
    end
end

begin # Plot horizontal velocity solutions versus simulation values
    plot(title="Layer Horizontal Velocities")
    plot!([data["layers"][k]["ys"][1,:] for k in layer_keys], [data["layers"][k]["vy"][1,:] for k in layer_keys], label=permutedims(["True $(k)" for k in layer_keys]), palette = :Dark2_5)

    
    ws = 0e3:100:20e3
    plot!(ws, [(res["solution"](ws) .* Vys_interp(ws)) for res in results], label=permutedims(["Solution for $(res["selected_layers"])" for res in results]), linestyle=:dash, palette = :Dark2_5)
    
    xlims!(0, 15e3)
    ylims!(0, 25)

    plot!(legend=:topleft)
end

begin # Plot cross-section view of results
    w_end = 15e3
    ws = 0e3:50:w_end

    zl = (400,1100)
    vl = (0, 20)

    # boundaries
    base_z = Spline1D(data["slices"]["center"]["ys"][1,:], data["slices"]["center"]["base"][:,1])
    surf_z = Spline1D(ys, data["layers"]["l0"]["geometry"][1,:])

    plot(ws/1e3, surf_z(ws), color="grey", linewidth=5, label="Surface")
    plot!(ws/1e3, base_z(ws), color="black", linewidth=5, label="Bed")
    
    velocity_end = zeros((length(results),))
    z_end = zeros((length(results),))

    for (idx, res) in enumerate(results)
        l1 = res["l1_z"](ws)
        l2 = res["l2_z"](ws)
        z = (l1 + l2) / 2 # z
        v = res["solution"](ws) .* Vys_interp(ws)
        scatter!(ws/1e3, z, zcolor=v, clims=vl, markerstrokewidth=0, linewidth=10, label="", vmin=0, vmax=15)

        z_end[idx] = z[end]
        velocity_end[idx] = v[end]
    end

    p1 = plot!(legend=:topright,
                colorbar_title="Horizontal Velocity [m/yr]",
                xlabel="Horizontal distance [km]", ylims=zl,
                ylabel="Elevation [m]")

    p2 = plot(velocity_end, z_end, label="Inversion result", markershape=:circle)
    plot!(xlabel="Horizontal Velocity [m/yr]", ylims=zl, xlims=vl)

    # Add "truth" from simulation results
    true_z = []
    true_vy = []
    for lk in layer_keys
        l = data["layers"][lk]
        idx = argmin(abs.(l["ys"][1,:].-w_end))
        push!(true_z, l["geometry"][1,idx])
        push!(true_vy, l["vy"][1,idx])
    end
    plot!(true_vy, true_z, label="True", linecolor=:gray, linestyle=:dash)

    plot(p1, p2, layout=grid(1, 2, widths=(0.7,0.3)), size=(800,(5/16)*800), dpi=600)

end
