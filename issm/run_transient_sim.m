function md = run_transient_sim(basal_friction_coefficient, rheology_n, rheology_paterson_C)
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

    md.timestepping=timesteppingadaptive();
    md.timestepping.time_step_max=20;
    md.timestepping.time_step_min=0.01;
    md.timestepping.start_time=0;
    md.timestepping.final_time=1000;
    
    md = solve(md, 'transient');
end