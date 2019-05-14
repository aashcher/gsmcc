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
%%
function [SV] = gsmcc1D(SVI, maxit, tol)
    global pol;

    V = calc_inc1D(SVI);
    W = gmres(@mul_A1D,V,[],tol,maxit);
%    W = bicgstab(@mul_A1D,V,tol,maxit);
    if strcmp(pol,'TE')
        VO = calc_radb1D(calc_emi1D_e(W));
    elseif strcmp(pol,'TM')
        VO = calc_radb1D(calc_emi1D_h(W));
    end
    SV = add_inc1D(SVI,VO);
end