function smb = smb_distribution(r)
    constant_radius = 20e3;
    ela_radius = 30e3;
    center_smb = 1;
    
    smb = zeros(size(r));
    r_inner = r(r < constant_radius);
    smb(r < constant_radius) = center_smb * ones(size(r_inner));
    r_linear = r(r >= constant_radius);
    smb(r >= constant_radius) = center_smb + ((r_linear - constant_radius) * (-1 * center_smb / (ela_radius - constant_radius)));
end