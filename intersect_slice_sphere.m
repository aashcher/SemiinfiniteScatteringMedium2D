%{
Copyright © 2022 Alexey A. Shcherbakov. All rights reserved.
This file is part of the SemiinfiniteScatteringMedium2D project.

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This sortware is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

mailto: alexey.shcherbakov@metalab.ifmo.ru
%}
%%
function [alp,pos] = intersect_slice_sphere(radii, cell_size, ncell_x, ...
																			nslice_per_cell, ind_cell, ind_slice)
    alp = zeros(ncell_x, 1);
    pos = zeros(ncell_x, 1);
    z = cell_size * ((nslice_per_cell - ind_slice + 0.5)/nslice_per_cell - 0.5);
    for i = 1 : ncell_x
        pos(i) = (i-0.5)/ncell_x - 0.5;
        R = radii(ind_cell, i);
        if abs(z) < R
            alp(i) = 2 * sqrt(R^2 - z^2) / (cell_size * ncell_x);
        end
    end
end