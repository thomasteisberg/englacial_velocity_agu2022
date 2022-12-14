clear all; close all;

%% Rheology, basal friction parameters

basal_friction_coefficient = 100;
rheology_n = 3.5;
rheology_paterson_C = -10;

%% Run ISSM simulation

md = run_transient_sim(basal_friction_coefficient, rheology_n, rheology_paterson_C);

%% Extract gridded slices from results

sol = md.results.TransientSolution(end);

y_spacing = 50; z_spacing = 50;
slice = extract_gridded_slice(md, sol, 0, y_spacing, z_spacing, 2000);
slice_offset_pos = extract_gridded_slice(md, sol, 50, y_spacing, z_spacing, 2000);
slice_offset_neg = extract_gridded_slice(md, sol, -50, y_spacing, z_spacing, 2000);

magnitude_grid_1 = sqrt(slice.Vz.^2 + slice.Vx.^2 + slice.Vy.^2);

figure
% Velocity magnitude color background
h = pcolor(slice.Y, slice.Z, magnitude_grid_1);
set(h, 'EdgeColor', 'none');
h.FaceColor = 'interp';
h.DisplayName = '';
axis equal
% Quiver plot of velocity directions
y_spacing = 20; z_spacing = 2;
hold on
quiver( slice.Y(1:y_spacing:end, 1:z_spacing:end), ...
        slice.Z(1:y_spacing:end, 1:z_spacing:end), ...
        slice.Vy(1:y_spacing:end, 1:z_spacing:end), ...
        slice.Vz(1:y_spacing:end, 1:z_spacing:end), 0.1, 'DisplayName', '');
% Surface and bed
plot(slice.ys, slice.surf, 'r', 'LineWidth', 2, 'DisplayName', '');
plot(slice.ys, slice.base, 'k', 'LineWidth', 2, 'DisplayName', '');

title('Cross-section velocity profile')
xlabel('Y [m]'); ylabel('Z [m]');
c = colorbar;
ylabel(c, 'Velocity Magnitude [m/yr]');
hold off
xlim([0, 40e3])
caxis([0, 30])

%% Synthetic layers
% Isochrones assuming end result is steady state
% Note that the technique does not depend upon any properties of how layers
% are physically formed, so not being completely physically correct layers
% shouldn't matter too much

mask = ~isnan(slice.Vx);

vx_interp = scatteredInterpolant(slice.Y(mask), slice.Z(mask), slice.Vx(mask));
vy_interp = scatteredInterpolant(slice.Y(mask), slice.Z(mask), slice.Vy(mask));
vz_interp = scatteredInterpolant(slice.Y(mask), slice.Z(mask), slice.Vz(mask));

vx_interp_offset_pos = scatteredInterpolant(slice_offset_pos.Y(mask), slice_offset_pos.Z(mask), slice_offset_pos.Vx(mask));
vy_interp_offset_pos = scatteredInterpolant(slice_offset_pos.Y(mask), slice_offset_pos.Z(mask), slice_offset_pos.Vy(mask));
vz_interp_offset_pos = scatteredInterpolant(slice_offset_pos.Y(mask), slice_offset_pos.Z(mask), slice_offset_pos.Vz(mask));

vx_interp_offset_neg = scatteredInterpolant(slice_offset_neg.Y(mask), slice_offset_neg.Z(mask), slice_offset_neg.Vx(mask));
vy_interp_offset_neg = scatteredInterpolant(slice_offset_neg.Y(mask), slice_offset_neg.Z(mask), slice_offset_neg.Vy(mask));
vz_interp_offset_neg = scatteredInterpolant(slice_offset_neg.Y(mask), slice_offset_neg.Z(mask), slice_offset_neg.Vz(mask));

particles_y = squeeze(slice.ys)';
particles_z = squeeze(slice.surf);

isochrone_years = [0 100 300 500 700 900 1100 1300];
timestep_years = 10;

layers = struct;

hold on

next_isochrone_index = 1;
for t = 0:timestep_years:max(isochrone_years)
    % Plot isochrone if requested
    if t >= isochrone_years(next_isochrone_index)
        mask = ~isnan(particles_z);
        layer_z = interp1(particles_y(mask), particles_z(mask), slice.ys);
        layer_name = ['l' num2str(t)];
        layers.(layer_name) = struct;
        layers.(layer_name).geometry = layer_z;
        layers.(layer_name).isochrone_year = t;
        
        mask = ~isnan(layer_z);
        layers.(layer_name).vx = vx_interp(slice.ys(mask), layer_z(mask));
        layers.(layer_name).vy = vy_interp(slice.ys(mask), layer_z(mask));
        layers.(layer_name).vz = vz_interp(slice.ys(mask), layer_z(mask));
        
        layers.(layer_name).ys = slice.ys(mask);
        layers.(layer_name).zs = layer_z(mask);
        
        layers.(layer_name).vx_offset_pos = vx_interp_offset_pos(slice.ys(mask), layer_z(mask));
        layers.(layer_name).vy_offset_pos = vy_interp_offset_pos(slice.ys(mask), layer_z(mask));
        layers.(layer_name).vz_offset_pos = vz_interp_offset_pos(slice.ys(mask), layer_z(mask));
        
        layers.(layer_name).vx_offset_neg = vx_interp_offset_neg(slice.ys(mask), layer_z(mask));
        layers.(layer_name).vy_offset_neg = vy_interp_offset_neg(slice.ys(mask), layer_z(mask));
        layers.(layer_name).vz_offset_neg = vz_interp_offset_neg(slice.ys(mask), layer_z(mask));
        
        plot(slice.ys, layer_z, 'DisplayName', ['Isochrone Yr ' num2str(t)])
        next_isochrone_index = next_isochrone_index + 1;
    end
    
    % Move particles
    delta_y = timestep_years * vy_interp(particles_y, particles_z);
    delta_z = timestep_years * vz_interp(particles_y, particles_z);
    particles_y = particles_y + delta_y;
    particles_z = particles_z + delta_z; 
end

hold off

%% Data Export

output_data = struct;
output_data.md = md;
output_data.layers = layers;
output_data.slices.center = slice;
output_data.slices.offset_pos = slice_offset_pos;
output_data.slices.offset_neg = slice_offset_neg;
filename = sprintf('sim_data_n%.1f_fric%d_Btemp%d.mat', rheology_n, basal_friction_coefficient, rheology_paterson_C)
save(filename, 'output_data', '-v7.3', '-nocompression') % Save in HDF5 format

