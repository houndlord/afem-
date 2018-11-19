function error4e = error4eP1Energy(c4n,n4e,gradExact,uApprox)
%% error4eP1Energy - compute the energy error on elements
%  Input:  c4n,n4e     - mesh
%          gradExact   - the exact gradient
%          uApprox     - AFEM P1 solution
%
%  Output: error4e     - the exact squared energy error on each element


    %% Compute grad U
    nrElems = size(n4e,1);
    gradU = zeros(nrElems,2);

    for elem = 1 : nrElems
        grads = [1,1,1;c4n(n4e(elem,:),:)']\[0,0;eye(2)];
        gradU(elem,:) = uApprox(n4e(elem,:))' * grads;
    end

    %% Compute the error
    
    error4e = integrate(c4n,n4e,@(n4p, Gpts4p, Gpts4ref) (...
	        sum((gradExact(Gpts4p) - gradU ).^2,2)),6);                    
    
    
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
