function slice = extract_gridded_slice(md, sol, x_slice, y_spacing, z_spacing, z_max)
%EXTRACT_GRIDDED_SLICE Extract 3D velocity vectors from a Y-Z slice of a
%TransientSolution
%   md - ISSM model object
%   sol - A TransientSolution output from ISSM
%   x_slice - X-location of the slice to extract [m]
%   y_spaicng, z_spacing - Resolution of gridded output in Y and Z [m]
%   z_max - Fixed maximum Z distance to include in grid. Set to -1 to
%   automatically select
arguments
    md
    sol {struct}
    x_slice double = 0
    y_spacing double = 100
    z_spacing double = 50
    z_max double = -1
end

% Interpolate results to grid on Y-Z cross section through the middle of the
% ice sheet
ys = min(md.mesh.y):y_spacing:max(md.mesh.y);

if z_max < 0
    z_max = max(sol.Surface)+100;
end

zs = 0:z_spacing:z_max;
[xq,yq,zq] = meshgrid(x_slice, ys, zs);

% Rescale Z points to match new thickness
rescaled_z = ((md.mesh.z - md.geometry.bed) .* (sol.Thickness ./ md.geometry.thickness)) + md.geometry.bed;

% Velocity field
Vz = squeeze(griddata(md.mesh.x, md.mesh.y, rescaled_z, sol.Vz, xq, yq, zq));
Vx = squeeze(griddata(md.mesh.x, md.mesh.y, rescaled_z, sol.Vx, xq, yq, zq));
Vy = squeeze(griddata(md.mesh.x, md.mesh.y, rescaled_z, sol.Vy, xq, yq, zq));
% Surface
surf_grid = squeeze(griddata(md.mesh.x, md.mesh.y, rescaled_z, sol.Surface, xq, yq, zq));
surf = median(surf_grid, 2, 'omitnan');
% Bed
base_grid = squeeze(griddata(md.mesh.x, md.mesh.y, rescaled_z, sol.Base, xq, yq, zq));
base = median(base_grid, 2, 'omitnan');

Vz(zs < base) = nan;
Vx(zs < base) = nan;

slice = struct;
slice.ys = squeeze(ys); slice.zs = squeeze(zs);
slice.Vx = Vx; slice.Vy = Vy; slice.Vz = Vz;
slice.surf = surf; slice.base = base;
slice.Y = squeeze(yq); slice.Z = squeeze(zq);
slice.x_slice = x_slice;

end

