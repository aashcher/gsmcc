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
function init_SI1D(es, eb, ec)
    global SI1;
    global SI2;
    global pol;
    
    SI1 = get_SI(es,eb,pol);
    SI2 = get_SI(eb,ec,pol);
%{    
    global no;
    global ky;
    global kz;

    SI1 = zeros(2,2,no);
    SI2 = zeros(2,2,no);
    
    if strcmp(pol,'TE')
        for i=1:no
            k1 = sqrt(es - ky(1,i)*ky(1,i));
            k2 = sqrt(ec - ky(1,i)*ky(1,i));
            SI1(1,1,i) = (k1 - kz(1,i))/(k1 + kz(1,i));
            SI2(1,1,i) = (kz(1,i) - k2)/(k2 + kz(1,i));
        end
    elseif strcmp(pol,'TM')
        for i=1:no
            k1 = sqrt(es - ky(1,i)*ky(1,i));
            k2 = sqrt(ec - ky(1,i)*ky(1,i));
            SI1(1,1,i) = (eb*k1 - es*kz(1,i))/(eb*k1 + es*kz(1,i));
            SI2(1,1,i) = (ec*kz(1,i) - eb*k2)/(eb*k2 + ec*kz(1,i));
        end
    else
        quit;
    end
    SI1(1,2,:) = 1 + SI1(1,1,:);
    SI1(2,2,:) = -SI1(1,1,:);
    SI1(2,1,:) = 1 + SI1(2,2,:);
    SI2(1,2,:) = 1 + SI2(1,1,:);
    SI2(2,2,:) = -SI2(1,1,:);
    SI2(2,1,:) = 1 + SI2(2,2,:);
%}
end

