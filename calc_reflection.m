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
% PM - power scattering matrix for all diffraction channels
% pa - power scattered to a given aperture for a given incidenct wave
function [pa, PM] = calc_reflection(structure, method, problem)
	pa = [];
	 % kz projections for the reciprocal lattice
	[kz] = fmm_kxz(method.no, problem.kx0, 0, problem.wavelength / structure.period, ...
									structure.eps_cov, structure.eps_cov);
	kzp = real(kz); % harmonics which propagate in the covering medium
	nop = nnz(kzr); % number of such harmonics
	indp = find(kzp); % indices of non-zero elements
	inda = acos(kzp) <= asin(problem.aperture); % inidices of harmonics, which constribute the given aperture
		% accumulate S-matrix for a given scattering structure:
	if strcmp(structure.type, 'cylinders')
			% initial scattering matrix:
		SM = calc_SMD_interface(method.no, problem.kx0, problem.wavelength / structure.period, ...
														structure.eps_med, structure.eps_cov, problem.polarization);
			% loop over vertical cells:
		for ih = 1 : structure.ncell_z
			for is = 1 : structure.nslice_per_cell % loop over slices in cell
					% parameters of piesewice-constant permittivity function
				[alp,pos] = intersect_slice_sphere(MR, structure.cell_size, ...
												structure.ncell_x, structure.nslice_per_cell, ih, is);
					% Fourier vector of the permittivity function:
				FE = calc_emn_bin(method.no, alp, pos, structure.eps_cyl(ih, :), structure.eps_med);
				SMt = fmm(method.no, problem.kx0, problem.wavelength / structure.period, ...
									structure.dh * problem.wavevector, ...
									FE, problem.polarization); % calculate slice S-matrix
				SM = mul_SM(SMt, SM);
			end
				% current power going to the given aperture:
			Vsca = smatrix_diffract(SM, problem.Vinc);
			pa = [pa; sum(abs(Vsca(inda, 2)^2), 'all')];
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif strcmp(structure.type, 'voronoi')
					% initial scattering matrix:
		SM = calc_SMD_interface(method.no, problem.kx0, problem.wavelength / structure.period, ...
														structure.eps_cov, structure.eps_cov, problem.polarization);
			% loop over slices:
		for is = 1 : structure.nslice % loop over slices in cell
			zs = structure.size_z * (1 - (is - 0.5) / structure.nslice);
				% parameters of piesewice-constant permittivity function
			[alp, pos, eps_sl] = intersect_slice_voronoi(zs, structure.vertices, ...
														structure.pvertices, structure.zlim, structure.eps);
				% Fourier vector of the permittivity function:
			FE = calc_emn_bin(method.no, alp, pos, eps_sl, 1);
			SMt = fmm(method.no, problem.kx0, problem.wavelength / structure.period, ...
								structure.dh * problem.wavevector, ...
								FE, problem.polarization); % calculate slice S-matrix
			SM = mul_SM(SMt, SM);
				% current power going to the given aperture:
			Vsca = smatrix_diffract(SM, problem.Vinc);
			pa = [pa; sum(abs(Vsca(inda, 2)^2), 'all')];
		end
	else
		error('calc_reflection: unknown structure type');
	end
		% power reflection matrix:
	PM = abs(SM(indp, indp, 2, 2)^2);
end