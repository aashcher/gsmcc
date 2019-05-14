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
function [VE] = calc_eff(VI, VD, eps1, eps2)
    global no;
    global pol;
    global ky;

    VE = zeros(2,no);
    Pi = 0;
      % accumulte incident and diffracted field power for each diffraction
      % order:
    if (strcmp(pol,'TE'))
        for i=1:no
            kz1 = sqrt(eps1 - ky(1,i)*ky(1,i));
            kz2 = sqrt(eps2 - ky(1,i)*ky(1,i));
            Pi = Pi + abs(VI(1,i)*VI(1,i))*real(kz1) + abs(VI(2,i)*VI(2,i))*real(kz2);
            VE(1,i) = abs(VD(1,i)*VD(1,i))*real(kz1);
            VE(2,i) = abs(VD(2,i)*VD(2,i))*real(kz2);
        end
    else
        for i=1:no
            kz1 = sqrt(eps1 - ky(1,i)*ky(1,i));
            kz2 = sqrt(eps2 - ky(1,i)*ky(1,i));
            Pi = Pi + abs(VI(1,i)*VI(1,i))*real(kz1/eps1) + abs(VI(2,i)*VI(2,i))*real(kz2/eps2);
            VE(1,i) = abs(VD(1,i)*VD(1,i))*real(kz1/eps1);
            VE(2,i) = abs(VD(2,i)*VD(2,i))*real(kz2/eps2);
        end
    end
    if (abs(Pi) > 1e-15)
        Pi = 1/Pi;
        for i=1:no
            VE(1,i) = VE(1,i)*Pi; VE(2,i) = VE(2,i)*Pi;
        end
    end
end

