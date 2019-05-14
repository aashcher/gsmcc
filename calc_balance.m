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
function [b] = calc_balance(VI, VD, eps1, eps2)
    global no;
    global pol;
    global ky;

    Pi = 0; Pd = 0;
    if (strcmp(pol,'TE'))
        for i=1:no
            kz1 = sqrt(eps1 - ky(1,i)*ky(1,i));
            kz2 = sqrt(eps2 - ky(1,i)*ky(1,i));
            Pi = Pi + abs(VI(1,i)*VI(1,i))*real(kz1) + abs(VI(2,i)*VI(2,i))*real(kz2);
            Pd = Pd + abs(VD(1,i)*VD(1,i))*real(kz1) + abs(VD(2,i)*VD(2,i))*real(kz2);
        end
    else
        for i=1:no
            kz1 = sqrt(eps1 - ky(1,i)*ky(1,i));
            kz2 = sqrt(eps2 - ky(1,i)*ky(1,i));
            Pi = Pi + abs(VI(1,i)*VI(1,i))*real(kz1/eps1) + abs(VI(2,i)*VI(2,i))*real(kz2/eps2);
            Pd = Pd + abs(VD(1,i)*VD(1,i))*real(kz1/eps1) + abs(VD(2,i)*VD(2,i))*real(kz2/eps2);
        end
    end
    if (abs(Pi) > 1e-15)
        b = abs(Pd/Pi-1);
    else
        b = Pd;
    end
end
