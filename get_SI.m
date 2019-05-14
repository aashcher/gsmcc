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
function S = get_SI(e1, e2, pol)
    global no;
    global ky;
    
    S = zeros(2,2,no);
    
    if strcmp(pol,'TE')
        for i=1:no
            k1 = sqrt(e1 - ky(1,i)*ky(1,i));
            k2 = sqrt(e2 - ky(1,i)*ky(1,i));
            S(1,1,i) = (k1 - k2)/(k1 + k2);
        end
    elseif strcmp(pol,'TM')
        for i=1:no
            k1 = sqrt(e1 - ky(1,i)*ky(1,i));
            k2 = sqrt(e2 - ky(1,i)*ky(1,i));
            S(1,1,i) = (e2*k1 - e1*k2)/(e2*k1 + e1*k2);
        end
    else
        disp('error in function get_SI: unknown polarization');
        quit;
    end
    
    S(1,2,:) = 1 + S(1,1,:);
    S(2,2,:) = -S(1,1,:);
    S(2,1,:) = 1 + S(2,2,:);
end
