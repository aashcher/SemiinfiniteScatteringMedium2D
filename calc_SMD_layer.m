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
%% function that calculates S-matrix of a homogeneous layer
%% input parameters:
% no: number of Fourier harmonics
% eps: layer medium permittivity
% kx0: zero harmonic wavevector projection
% kg: wavevector step
% kh: layer thickness multiplied by the wavenumber in vacuum
%% output parameters:
% SML - S-matrix
%% implementation:
function [SMD] = calc_SMD_layer(no, kx0, kg, kh, eps)
	kz = fmm_kxz(no, kx0, 0, kg, eps, eps);
	SMD = zeros(no,2,2);
	SMD(:,1,2) = exp((1i*kh)*kz);
	SMD(:,2,1) = SMD(:,1,2);
end