function [L2norm4Omega, L2norm4p] = L2Norm(c4n,n4e,f,degree)
%% L2Norm - compute the L2 norm of a function on a triangulation
%  Input:  c4n,n4e     - mesh
%          f           - function, the norm of which is computed;
%                        has to suit input/output behaviour needed by integrate
%          degree      - degree of accuracy used in the quadrature formula
%  Output: L2norm4p    - ||f||_(L^2(Part))
%          L2normOmega - ||f||_(L^2(Omega))

    %% Compute the L2 norm.
    L2norm4p = integrate(c4n, n4e, ...
    		@(n4p,Gpts4p,Gpts4ref)(f(n4p,Gpts4p,Gpts4ref)).^2, degree);
    L2norm4Omega = sqrt(sum(L2norm4p));
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
