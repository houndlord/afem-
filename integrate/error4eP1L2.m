function error4e = error4eP1L2(c4n,n4e,uExact,uApprox)
%% error4eP1L2 - compute the exact error on elements
%  Input:  c4n,n4e     - mesh
%          uExact      - the exact function
%          uApprox     - AFEM P1 solution
%
%  Output: error4e     - the exact squared error on each element

%% Compute error
  error4e = (integrate(c4n,n4e,@(n4p, Gpts4p, Gpts4ref) (...
	        uExact(Gpts4p)-...
            ((1 - Gpts4ref(:,1) - Gpts4ref(:,2))*uApprox(n4p(:,1)) +...
            (Gpts4ref(:,1))*uApprox(n4p(:,2)) +...
            (Gpts4ref(:,2))*uApprox(n4p(:,3)))).^2,4));      
        
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
