function [eta4s,n4s] = estimateRT0EtaSides(f, g, u4Db, p, u, c4n, n4e, n4sDb, n4sNb)
%% estimateRT0EtaSides - error estimator for RT0 element
% Estimate the energy error of the RT0 finite element solution by the
% jumps of the discrete solution's P0 component along the sides.
%
% Input:     f	       right-hand side of the problem definition
%            g         Neumann boundary condition
%            u4Db      Dirichlet boundary condition
%            x         RT0 basis coefficients of (u,p) given by solve
%            c4n       coordinates for the nodes of the mesh
%            n4e       nodes for the elements of the mesh
%            n4sDb	   the nodes of the sides in the Dirichlet boundary
%            n4sNb	   the nodes of the sides in the Neumann boundary
%
% Output:    eta4s     error for each side
%            n4s       nodes of the sides for which the error was computed

    %% initialisation
    s4n = computeS4n(n4e);
    n4s = computeN4s(n4e);
    e4s = computeE4s(n4e);

    length4s = computeLength4s(c4n,n4s);
    area4e   = computeArea4e(c4n,n4e);

    if size(n4sDb,1)>0
        s4Db = diag(s4n(n4sDb(:,1),n4sDb(:,2)));
    else
        s4Db = [];
    end

    s4Nb = diag(s4n(n4sNb(:,1),n4sNb(:,2)));

    %% jump and volume term
    jump4s = integrate(c4n, n4s, @(parts, Gpts4p, Gpt4ref)...
         computeJump4s(parts, Gpts4p, Gpt4ref, p, n4e, c4n,...
         e4s, s4Db, s4Nb, n4s, length4s, u4Db), 2);

    vol4e = integrate(c4n, n4e, @(parts, Gpts4p, Gpt4ref) f(Gpts4p).^2, 2);

    TPlus4s  = e4s(:,1);
    TMinus4s = e4s(:,2);
    vol4eTMinus = zeros(size(n4s,1),1);
    vol4eTMinus(TMinus4s(TMinus4s~=0)) = vol4e(TMinus4s(TMinus4s~=0));
    area4eTMinus = zeros(size(n4s,1),1);
    area4eTMinus(TMinus4s(TMinus4s~=0)) = area4e(TMinus4s(TMinus4s~=0));

    eta4s = length4s .* jump4s + area4e(TPlus4s).*vol4e(TPlus4s)/3 ...
            + area4eTMinus.*vol4eTMinus/3;

end

%% Jumps
function jumps4s = computeJump4s(n4parts, Gpts4p, Gpt4ref, p, n4e, c4n, ...
		                 e4s,s4Db,s4Nb, n4s, length4s, u4Db)
    tangent4s = computeTangent4s(c4n,n4parts); 
    mid4e = computeMid4e(c4n, n4e);
    
    % components of the solution
    p1 = p(:,1);
    p2 = p(:,2);
    p3 = p(:,3);
    
    TPlus4s = e4s(:,1);
    TMinus4s = e4s(:,2);
    %Use element 1 for calculation and overwrite that later with the correct value
    TMinus4s(TMinus4s == 0) = 1; 

    valTPlus =  sum ( (p1(TPlus4s) * [1 0]  + p2(TPlus4s) * [0 1] ...
                       + (Gpts4p - mid4e(TPlus4s,:)) ...
		       .* [p3(TPlus4s) p3(TPlus4s)]).* tangent4s, 2);
    valTMinus = sum ( (p1(TMinus4s) * [1 0] + p2(TMinus4s) * [0 1] ...
                       + (Gpts4p - mid4e(TMinus4s,:)) ...
		       .* [p3(TMinus4s) p3(TMinus4s)]).* tangent4s, 2);
    % p' in t
    valTMinus(s4Db) = (u4Db(c4n(n4s(s4Db,2),:)) - u4Db(c4n(n4s(s4Db,1),:))) ...
                    ./ length4s(s4Db);
    jumps4s = (valTPlus - valTMinus) .^ 2;
    jumps4s(s4Nb) = 0;
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
