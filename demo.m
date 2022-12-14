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
clc;
format long;
%%
%{
The code calculates 2D scattering of a polarized electromangetic wave by a
semi-infinite structure, which is represented either as a set of random
cylinders or a Voronoi diagram with each cell having some prescribed
dielectric permittivity.
             ___          ^ z
\	\           / |         |
 \ \         / /|         |
  \ \|      / /           |     # incident and reflected beam
	___|     / /            |
                          |
==============================> x
o O O   o o  o O o O   O O o
 o   O o   o o O O   oo   O O   # semi-infinite scattering medium
O   O oo OO   o   o o O O o     # with wavelength-scale inhomogneities


The supercell approach is used, which means that the medium is simulated as
being periodic in the horizontal X direction. A numerical limit for reflection
characteristics for increasing period can be calculated. Enseble averaged
power reflection can be calculated.

The periodic problem is solved by the Fourier Modal Method. For this pupose
a given realization of scattering structure (a sufficiently thick layer 
-H < z < 0 filled with inclusions of different permittivities) is divided 
into a set of plane slices, which thickness dh is small comparable to the
wavelength, and the spatial permittivity distribution in each slice is
represented as a piecewise-constant function of the coordinate X in a
period. FMM is used to find the S-matrix of each slice and the S-matrix
algorithm is used to calculate an S-matrix of a layer of N_s slices.
Convergence of the power reflection in the dependence of N_s can be traced.

In case of random cylinders to generate a structure a thick period is
divided into square 'cells', each cell can either contain a cylinder of
random radius or be empty. This is an example and any other structure 
generation algorithm preserving the periodicity can be used.

In case of Voronoi structure

For reliable simulations the number of Fourier harmonics used in the FMM
should be at least about four time larger than the period-to-wavelength
ratio.
%}

	%% initialization of a scattering structure
	% in case of random spheres:
structure.type = 'cylinders';
structure.cell_size = 0.5; % size of one "cell" / dd
structure.nslice_per_cell = 20; % number of slices to divide a cell / ndd
structure.ncell_x = 15; % number of cells in horizontal dimension (in one period) / nd_p
structure.ncell_z = 50; % number of cells in vertical dimension / nd_h
structure.dh = structure.cell_size / structure.nslice_per_cell; % slice thickness / dh
 % minimum and maximum radii for cylinders
structure.rmin = 0.3 * structure.cell_size;
structure.rmax = 0.5 * structure.cell_size; % cylinders cannot intersect
 % array of booleans determining whether a cell is empty or not:
structure.radii = rand(structure.ncell_z, structure.ncell_x) > .4;
 % randomly generate radii:
structure.radii = structure.radii .* (rand(structure.ncell_z, structure.ncell_x) * ...
											(structure.rmax - structure.rmin) + structure.rmin);
 % period of the supercell:
structure.period = structure.cell_size * structure.ncell_x;
 % material parameters:
structure.eps_cyl = (4+0.01*1i) * ones(structure.ncell_z, structure.ncell_x); % permittivities of cylinders
structure.eps_med = 2.2; % medium, which surrounds the cylinders
structure.eps_cov = 1; % permittivity of semi-infinite medium in contact with the structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
	% in case of Voroni diagram medium:
structure.type = 'voronoi';
structure.size_x = 5; % 
structure.size_z = 20; % 
structure.period = structure.size_x; % 
structure.nslice = 400; % number of slices
structure.dh = structure.size_z / structure.nslice;
	% generate Voronoi diagram:
num_points = 100; % number of randomly scattered points in one period
points_X = structure.size_x * rand(num_points, 1);
points_Z = structure.size_z * rand(num_points, 1);
 % generate points for three periods to make the diagram periodic:
points_X = [points_X; (points_X-structure.size_x); (points_X+structure.size_x)];
points_Z = [points_Z; points_Z; points_Z];
[V, C, XY] = VoronoiLimit(points_X, points_Z, 'figure', 'off');
structure.vertices = V;
structure.pvertices = C;
 % determine maximum and minimum z-coordinates for each polygon:
structure.zlim = zeros(length(C), 2);
for i = 1 : length(C)
	structure.zlim(i,1) = min(V(C{i,1}(:),2));
	structure.zlim(i,2) = max(V(C{i,1}(:),2));
end
 % randomly generate permittivities for each Voronoi cell for a given set of media:
structure.eps_media = [2.24+0.01*1i, 1.67, 1.95+0.02*1i]; % define permittivities of different materials within the period
structure.num_media = numel(structure.eps_media); % number of different media in a period
eps_vertmed = floor(rand(num_points, 1) * structure.num_media);
structure.eps = zeros(length(C), 1);
for im = 1 : structure.num_media
	ind = (eps_vertmed - i + 1) < 1e-12;
	structure.eps(ind) = structure.eps_media(i);
	structure.eps(ind + num_points) = structure.eps(ind);
	structure.eps(ind + 2*num_points) = structure.eps(ind);
end
structure.eps_cov = 1; % permittivity of semi-infinite medium in contact with the structure
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%% problem and method:
 % define parameters of radiation and output
problem.wavelength = 1; % wavelength
problem.wavevector = 2 * pi / problem.wavelength;
problem.angle = 1e-3; % angle of incidence in degrees
problem.kx0 = real(sqrt(structure.eps_cov)) * sin(problem.angle * pi / 180); % Bloch wavevector
problem.aperture = 0.3; % numerical aperture of receiving instrument
problem.polarization = 'TE'; % TE or TM

 % define parameters for the Fourier Modal Method:
method.no = 5*ceil(structure.period / problem.wavelength); % number of Fourier harmonics
method.ic = ceil(method.no / 2); % index of the central '0'-th harmonic

	% vector of incident wave amplitudes:
problem.Vinc = zeros(method.no,2);
problem.Vinc(method.ic, 2) = 1; % plane wave incidence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reflection calculation
% pa - array of  calculated power scattered to a given aperture with
% increaing structure thickness from 0 to maximum depth
% PM - power scattering matrix for all diffraction channels
[pa, PM] = calc_reflection(structure, method, problem);
