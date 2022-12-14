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
function [alp, xc, eps_sl] = intersect_slice_voronoi(zs, V, C, zlim, eps)
	alp = []; 
	xc = []; 
	eps_sl = [];
	for i=1:length(C)
		if (zs < zlim(i,2)) && (zs > zlim(i,1))
			xit = zeros(2,1);
			m = 1;
			nv = length(C{i,1});
			x2 = V(C{i,1}(nv,1),1);
			z2 = V(C{i,1}(nv,1),2);
			for j = 1 : nv
				x1 = x2; z1 = z2;
				x2 = V(C{i,1}(j,1),1);
				z2 = V(C{i,1}(j,1),2);
				if abs(z2 - z1) < 1e-15
					if abs(zs - z1) < 1e-15
						xit(1,1) = x1;
						xit(2,1) = x2;
						m = 3;
					end      
				elseif (zs <= max([z1,z2])) && (zs >= min([z1,z2]))          
					if m < 3
						xit(m,1) = x1 + (x2-x1)*(zs-z1)/(z2-z1);
						m = m + 1;
					end	
				end	
			end
			alp = [alp, abs(xit(2,1)-xit(1,1))];
			xc = [xc, 0.5*(xit(2,1)+xit(1,1))];
			eps_sl = [eps_sl, eps(i)];
		end
	end
end