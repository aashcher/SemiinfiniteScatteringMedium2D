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
%% description:
% Fourier matrix of the permittivity and the inverse permittivity of a
% 1D binary grating
%% input:
% no - number of Fourier harmonics
% alps: relative widths of the grating period elements (should meet the conditions
%   0 < alps(i) < 1 and sum(alps) == 1)
% poss: relative positions of the centers of the grating period elements
%   (should not overlap)
% eps: permittivities of the grating period elements
%% output:
% FE: Fourier matrix of size (2*no-1,2), which contains
% Fourier components of the permittivity (FE(:,1)) and the inverse
% permittivity (FE(:,2))
%% implementation
function [FE] = calc_emn_bin(no, alps, poss, eps, epsm)
	if (length(alps) ~= length(poss)) || (length(alps) ~= length(eps)) || (abs(sum(alps)) > 1)
		error("calc_emn_bin input error");
	end

	FE = zeros(2*no-1,2);
	ind = transpose(linspace(1,no-1,no-1));

	for k = 1:length(alps)
		ifun = sin(ind*pi*alps(k))./(ind*pi);
		te = exp(-(2*pi*1i*poss(k))*ind);
			% zero harmonics:
		FE(no,1) = FE(no,1) + eps(k)*alps(k);
		FE(no,2) = FE(no,2) + alps(k)/eps(k);
			% non-zero harmonics:
		tmp = eps(k)*ifun;
		FE(no+1:2*no-1,1) = FE(no+1:2*no-1,1) + tmp.*te;
		FE(no-1:-1:1,1) = FE(no-1:-1:1,1) + tmp.*conj(te);
		tmp = (1/eps(k))*ifun;
		FE(no+1:2*no-1,2) = FE(no+1:2*no-1,2) + tmp.*te;
		FE(no-1:-1:1,2) = FE(no-1:-1:1,2) + tmp.*conj(te);
	end
	wd = sum(alps, 'all');
  FE(no, 1) = FE(no, 1) + epsm * (1 - wd);
  FE(no, 2) = FE(no, 2) + (1 - wd) / epsm;
end
%
% end of calc_emn_bin
%