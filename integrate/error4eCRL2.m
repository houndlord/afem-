function error4e = error4eCRL2(c4n,n4e,uExact,uApprox)

s4p = computeS4e(n4e);

error4e = (integrate(c4n,n4e, ...
    @(n4p, Gpts4p, Gpts4ref) ((uExact(Gpts4p)-...
        ((-1 + 2*Gpts4ref(1) + 2*Gpts4ref(2))*uApprox(s4p(:,2)) +...
        (1-2*Gpts4ref(1))*uApprox(s4p(:,3)) +...
        (1-2*Gpts4ref(2))*uApprox(s4p(:,1)))).^2), 12));
end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009-2015
% Numerical Analysis Group
% Prof. Dr. Carsten Carstensen
% Humboldt-University
% Departement of Mathematics
% 10099 Berlin
% Germany
%
% This file is part of AFEM.
%
% AFEM is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% AFEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
