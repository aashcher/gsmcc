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
function calc_kyz(ky0, kg, eb)
    global no;
    global kz;
    global ky;
    
    ky = zeros(1,no);
    kz = zeros(1,no);
    ii = floor(no/2)+1;
    for i=1:no
        ky(1,i) = ky0 + kg*double(i-ii);
        kz(1,i) = sqrt(eb - ky(1,i)*ky(1,i));
%        if arg(kz(1,i)) < -1e-8
%            kz(1,i) = -kz(1,i);
%        end
    end
end