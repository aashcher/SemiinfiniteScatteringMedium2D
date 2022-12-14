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
%% function that calculates S-matrix of an interface between two homogeneous isotropic media
%% input parameters:
% no: number of Fourier harmonics
% eps1, eps2: medium permittivities below and above
% kx0: zero harmonic wavevector projection
% kg: wavevector step
% pol: polarization, either 'TE' or 'TM'
%% output parameters:
% SMD: diagonal interface S-matrix
%%
function [SMD] = calc_SMD_interface(no, kx0, kg, eps1, eps2, pol)
	[kz1, kz2] = fmm_kxz(no, kx0, 0, kg, eps1, eps2);

	SMD = zeros(no,2,2);
	if strcmp(pol,'TE')
		SMD(:,1,1) = (kz1-kz2)./(kz1+kz2);
		SMD(:,2,1) = 1 + SMD(:,1,1);
		SMD(:,2,2) = -SMD(:,1,1);
		SMD(:,1,2) = 1 + SMD(:,2,2);
	elseif strcmp(pol,'TM')
		SMD(:,1,1) = (eps2*kz1-eps1*kz2)./(eps2*kz1+eps1*kz2);
		SMD(:,2,1) = 1 + SMD(:,1,1);
		SMD(:,2,2) = -SMD(:,1,1);
		SMD(:,1,2) = 1 + SMD(:,2,2);
	else
		error('function calc_SM_interface: unknown polarization');
	end
end


