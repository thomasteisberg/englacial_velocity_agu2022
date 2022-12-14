
basal_friction_coefficient = 500;
rheology_n = 2.5;
rheology_paterson_C = -20;

%%

md = model;

% Mesh

md = triangle(md, 'domain.exp', 4000);

md = setmask(md, '', ''); % all grounded
md = parameterize(md, 'syn_glacier_constant_slope.par');
md = extrude(md, 10, 1.0);
md = setflowequation(md, 'HO', 'all');

% Update with fields being tested
md.friction.coefficient = basal_friction_coefficient * ones(md.mesh.numberofvertices,1);
md.materials.rheology_B = paterson((273 + rheology_paterson_C) * ones(md.mesh.numberofvertices,1));
md.materials.rheology_n = rheology_n * ones(md.mesh.numberofelements,1);

md.transient.isthermal = 0;

% md.timestepping.final_time = 100;
% md.timestepping.time_step = 1.0;
md.timestepping=timesteppingadaptive();
md.timestepping.time_step_max=20;
md.timestepping.time_step_min=0.01;
md.timestepping.start_time=0;
md.timestepping.final_time=500; % TODO: 1500

md = solve(md, 'transient');

%%

%% Plot evolution of ice thickness to try to visually check if we're in steady state

for i = 1:length(md.results.TransientSolution)
    plotmodel(md,'data',md.results.TransientSolution(i).Thickness);
    c = colorbar;
    ylabel(c, 'Thickness [m]');
    caxis([0 1500])
    view(2)
    title(['Timestep ' num2str(i)]); xlabel('X [m]'); ylabel('Y [m]');
    drawnow
end


%%

sol = md.results.TransientSolution(end);

plotmodel(md, 'data', sqrt(sol.Vx.^2 + sol.Vy.^2), 'caxis#1', [0 30], ...
    'data', sol.Thickness, ...
    'view#all', [0 90], 'xlabel#all', 'X', 'ylabel#all', 'Y' ...
    )

%%

y_spacing = 50;
surface_slice = extract_surface_velocity_line(md, sol, 0, y_spacing);

magnitude_surf = abs(surface_slice.Vy);

plot(surface_slice.ys, magnitude_surf)
xlim([0, 40e3])
