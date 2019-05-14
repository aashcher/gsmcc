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
function S = mul_SSD(S1, S2)
    TS = 1 - S1(2,2,:).*S2(1,1,:);
    TS = 1./TS;
    S(1,1,:) = S1(1,1,:) + S1(1,2,:).*S2(1,1,:).*S1(2,1,:).*TS;
    S(1,2,:) = S1(1,2,:).*S2(1,2,:).*TS;
    S(2,1,:) = S1(2,1,:).*S2(2,1,:).*TS;
    S(2,2,:) = S2(2,2,:) + S2(2,1,:).*S1(2,2,:).*S2(1,2,:).*TS;
end