%{
Copyright © 2019 Alexey A. Shcherbakov. All rights reserved.

This file is part of gsmcc.

gsmcc is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

gsmcc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sphereml. If not, see <https://www.gnu.org/licenses/>.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
gsmcc - the Generalized Source Method in Curvilinear Coordinates

This is the main file for the set of functions implementing the GSMCC
method for 1D gratings in the collinear diffraction case. 
Gratings can be sandwitched inside any planar multilayer
structure.

This code is intended for non-commercial use only.

If any output of the gsmcc or of a modified vesion of the current code is
published, please, refer to:
[1] A. A. Shcherbakov and A. V. Tishchenko, Efficient curvilinear coordinate
method for grating diffraction simulation, Opt. Express 21, 25236-25247 (2013)
[2] A. A. Shcherbakov, Curvilinear coordinate Generalized Source Method for 
gratings having continuous piecewise-smooth profiles arXiv:1904.05757
All the detailes regarding the method can be found in these publications.

If you have any questions, need a more efficient, flexible or parallel implementation,
please, write to alex.shcherbakov@gmail.com
%}
%%
clc;
format long;
%% GLOBAL VARIABLES / DO NOT MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global no; % number of Fourier harmonics in truncated series; zero order harmonics have index floor(no/2)+1
global ns; % number of slices for integration
global hh; % total calculation layer depth
global eb; % basis layer permittivity
global cwz; % matrix of slice coordinates and integration weigths
global ky; % vector of in-plane grating vector projections
global kz; % vector vertical wavevector prpojections in the basis medium
global SI1; % diagonal S-matrix of the part of the structure below the grating
global SI2; % diagonal S-matrix of the part of the structure above the grating
global pol; % polarization, either 'TE' or 'TM'
global MG;

%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pol = 'TM'; % polarization
gp = 0.5; % grating period
gh = 0.1; % grating depth (has no sense for multilayer gratings)
wl = 0.6; % wavelength
es = 2.25; % substrate permittivity
eg = 2.25 + 0*1i; % grating permittivity (has no sense for multilayer gratings)
ec = 1.0; % cover permittivity
wv = 2*pi/wl; % vacuum wavenumber
kh = wv*gh;
kg = wl/gp;
ab = 1.0; % a/b ratio, see the publications, better leave to be 1
th = 10.0; % incidence angle (in the Air)
ky0 = sin(pi*th/180); % incident wavevector projection (Bloch wavevector modulo)

eb = 0.5*(eg + ec); % basis grating layer permittivity (may affect the convergence)
no = 16; % number of Fourier harmonics in truncated series; zero order harmonics have index floor(no/2)+1
ns = 100; % number of slices for integration
hh = kh/ab; % grating region depth

%% INITIALIZATIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate wavevector projections for plane harmonics:
calc_kyz(ky0, kg, eb);
% calculate diagonal scattering matrices of structure under and above the grating:
EL1 = []; % permittivities of layers below the grating (leave empty vector if there no layers)
khL1 = []; % widths of layers below the grating normalized by the vacuum wavenumber
EL2 = []; % permittivities of layers above the grating
khL2 = []; % widths of layers above the grating normalized by the vacuum wavenumber
init_SSI1D(es, EL1, khL1, eb, EL2, khL2, ec);

% incident wave amplitude matrix:
% VI(1,:) - amplitudes of plane waves coming from below the structure
% VI(2,:) - amplitudes of plane waves coming from above the structure
VI = zeros(2,no);
VI(2,floor(no/2)+1) = 1; % zero order plane wave with unit amplitude

%% GRATINGS (uncomment one of the following options) %%%%%%%%%%%%%%%%%%%%%%
% *1* if simulating a simple sinusoidal grating:
init_sin1D(kg, ab, eg, ec);

% *2* if simulating a multilayer sinusoidal grating:
%{
ni = 2; % number of sinusoidal interfaces
ka = [kh*0.5, kh*0.5]; % amplitudes of interfaces normalized by the vacuum wavenumber (ni)
kt = [wv*0.05]; % vertical distance between zero-levels of the sinusoidal
% functions normalized by the vacuum wavenumber (ni-1)
% make sure that the interface functions do not intersect !
esl = [eg]; % permittivities of layers between interfaces (ni-1)
init_sinml1D(kg, ab, ni, ka, kt, esl, es, ec);
%}
% *3* if simulating a polygonal shape grating (example of trapezoidal profile):
%{
np = 6; % 6 points of the polyline interface
YZ = zeros(2,np);
% define y and z coordinates of vertices (from left right):
% y coordinates are in fractions of the period, from -0.5 to 0.5:
YZ(1,1) = -0.5; YZ(1,6) = 0.5; % the first and the last points should be exactly in the beginning and at the end of the pertiod
YZ(1,2) = -0.3; YZ(1,3) = -0.1; YZ(1,4) = 0.1; YZ(1,5) = 0.3;
% z coordinates are normalized to the vacuum wavenumber:
YZ(2,1) = -0.5*kh; YZ(2,2) = -0.5*kh; YZ(2,3) = 0.5*kh; YZ(2,4) = 0.5*kh; YZ(2,5) = -0.5*kh; YZ(2,6) = -0.5*kh;
% relative values of z points is unimportant, grating depth is
% calculated as a difference between maximum and minimum values
init_polyl1D(kg, YZ, ab, eg, ec);
%}
%% EXAMPLE OF A SINGLE RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itermax = 500; % maximum number of iterations for the GMRes
tol = 1.e-5; % convergence tolerance for the GMRes

% calculate diffraction amplitude vector:
% use the following function for sinusoidal shape gratings:
SV = gsmcc1D(VI,itermax,tol);
% use the following function for polygonal shape gratings:
%SV = gsmcc1Dc(VI,itermax,tol);
% output:
% SV(1,:) - amplitudes of diffracted harmonics below the structure in the substrate
% SV(2,:) - amplitudes of diffracted harmonics above the structure in the cover
SE = calc_eff(VI,SV,es,ec); % calculate corresponding matrix of diffraction efficiencies

disp(calc_balance(VI,SV,es,ec)); % power balance; should be near-zero for dielectric structures
disp(SE(1,floor(no/2)+1)); % zero order efficiency in the substrate
disp(SE(2,floor(no/2)+1)); % zero order efficiency in the cover

%% CONVERGENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% {
% for any new grating structure and parapmeter range one should check the convergence to be sure
% in the validity of results; value of "ns" should be large enough, and "tol"
% should be small enough for a good convergence
ns = 200;
no1 = 5; % starting numner of Fourier orders
no2 = 20; % end number of Fourier orders
nno = no2 - no1;
tol = 1e-10;
dataCONV = zeros(4,nno); % array to save the data
SV2 = zeros(2,no1);
for no = no1:no2
      % all initializations should be done for each new "no" since the
      % sizes of all vectors and matrices change
    calc_kyz(ky0, kg, eb);
    init_SSI1D(es, EL1, khL1, eb, EL2, khL2, ec);
    VI = zeros(2,no);
    VI(2,floor(no/2)+1) = 1; % zero order plane wave with unit amplitude
    init_sin1D(kg, ab, eg, ec);
    
    SV1 = SV2;
    SV2 = gsmcc1D(VI,itermax,tol);
    b = calc_balance(VI,SV2,es,ec);
    SE = calc_eff(VI,SV2,es,ec);
    
    if no > no1
          % difference between diffracted amplitudes of 0th harmonic:
        dVD = SV2(1,floor(no/2)+1) - SV1(1,floor((no-1)/2)+1);
        dataCONV(1,no-no1) = 1/no;
        dataCONV(2,no-no1) = norm(dVD);
        dataCONV(3,no-no1) = SE(1,floor(no/2)+1);
        dataCONV(4,no-no1) = b;
    end
end

plot(dataCONV(1,:),dataCONV(2,:));
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');   
%}
%% SPECTRAL CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
VI = zeros(2,no);
VI(2,floor(no/2)+1) = 1;

n_wl = 350; % number of wavengths
wl1 = 0.4; % initial wavelength
wl2 = 0.75; % final wavelenth
dwl = (wl2 - wl1)/n_wl; % wavelength step
dataSPEC = zeros(2,n_wl);

for iwl = 1:n_wl
    wl = wl1 + iwl*dwl; % current wavelength
    wv = 2*pi/wl; kg = wl/gp; kh = wv*gh;
    % update permittivities here if the dispersion is significant
    % ...
    
    calc_kyz(ky0, kg, eb);
    init_SSI1D(es, EL1, khL1, eb, EL2, khL2, ec); % if khL1,khL2 are not empty, they should be updated
    init_sin1D(kg, ab, eg, ec);
    SV = gsmcc1D(VI,itermax,tol);
    b = calc_balance(VI,SV,es,ec);
    SE = calc_eff(VI,SV,es,ec);
    
    dataSPEC(1,iwl) = wl;
    dataSPEC(2,iwl) = SE(1,floor(no/2)+1);
end

plot(dataSPEC(1,:),dataSPEC(2,:));
%}
