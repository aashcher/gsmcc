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
function [W] = mul_A1D(V)
    global pol;
    
    if strcmp(pol,'TE')
        W = V - calc_rad1D(calc_emi1D_e(V));
    elseif strcmp(pol,'TM')
        W = V - calc_rad1D(calc_emi1D_h(V));
    end
end