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
function init_SSI1D(es, EL1, khL1, eb, EL2, khL2, ec)
    global SI1;
    global SI2;
    global pol;
    
    if (length(EL1) ~= length(khL1)) || (length(EL2) ~= length(khL2))
        disp('error in function init_SSI1D: lengths of permittivity and width vectors are different');
        quit;
    end
    
    if isempty(EL1) == 1
        SI1 = get_SI(es,eb,pol);
    else
        nl = length(EL1);
        SI1 = get_SI(es,EL1(1,1),pol);
        SI1 = mul_SSD(SI1,get_SL(EL1(1,1),khL1(1,1)));
        for k = 1:(nl-1)
            SI1 = mul_SSD(SI1,get_SI(EL1(1,k),EL1(1,k+1),pol));
            SI1 = mul_SSD(SI1,get_SL(EL1(1,k+1),khL1(1,k+1)));
        end
        SI1 = mul_SSD(SI1,get_SI(EL1(1,nl),eb,pol));
    end

    if isempty(EL2) == 1
        SI2 = get_SI(eb,ec,pol);
    else
        nl = length(EL2);
        SI2 = get_SI(eb,EL2(1,1),pol);
        SI2 = mul_SSD(SI2,get_SL(EL2(1,1),khL2(1,1)));
        for k = 1:(nl-1)
            SI2 = mul_SSD(SI2,get_SI(EL2(1,k),EL2(1,k+1),pol));
            SI2 = mul_SSD(SI2,get_SL(EL2(1,k+1),khL2(1,k+1)));
        end
        SI2 = mul_SSD(SI2,get_SI(EL2(1,nl),ec,pol));
    end
end

