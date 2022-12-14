function b = base_elevation(r) 
    bed_elevation_top = 500; % m
    width_scale = 5e3;
    transition_to_linear = 500; % m
    
    linear_slope = -1*(4*exp(transition_to_linear/width_scale)*(-1 + exp(transition_to_linear/width_scale))*bed_elevation_top)/((1 + exp(transition_to_linear/width_scale))^3 * width_scale)
    intersect = 4*bed_elevation_top*exp(transition_to_linear/width_scale) ./ (1+exp(transition_to_linear/width_scale)).^2;
    intersect = intersect - linear_slope*transition_to_linear;
    
    b = zeros(size(r));
    r_smooth = r(r < transition_to_linear);
    b(r < transition_to_linear) = 4*bed_elevation_top*exp(r_smooth/width_scale) ./ (1+exp(r_smooth/width_scale)).^2;
    r_linear = r(r >= transition_to_linear);
    b(r >= transition_to_linear) = r_linear * linear_slope + intersect;
end