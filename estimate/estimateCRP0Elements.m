function eta4e = estimateCRP0Elements(f,g,u4Db,u,p,c4n,n4e,n4sDb,n4sNb,mu)

    if nargin < 10
        mu = 1;
    end
    %% initialisation
    n4s = computeN4s(n4e);
    s4e = computeS4e(n4e);

    nrElems = size(n4e,1);

    area4e = computeArea4e(c4n,n4e);
    mid4e = computeMid4e(c4n,n4e);

    length4s = computeLength4s(c4n,n4s);

    %% Compute the CR-gradient of u
    grad4e1 = zeros(nrElems,2);
    grad4e2 = zeros(nrElems,2);
    for elem = 1 : nrElems
        grads = [1,1,1;c4n(n4e(elem,:),:)'] \ [0,0;eye(2)];
        nc_grads = [-1,1,1;1,-1,1;1,1,-1] * grads;
        grad4e1(elem,:) = u(s4e(elem,:),1)'*nc_grads;
        grad4e2(elem,:) = u(s4e(elem,:),2)'*nc_grads;
    end
    
    %% jumps
%     jumpN1 = mu*P0NormalJump(c4n,n4e,n4sDb,n4sNb,grad4e1-[p,p],@(x)gNb(x,1,g));
%     jumpN2 = mu*P0NormalJump(c4n,n4e,n4sDb,n4sNb,grad4e2-[p,p],@(x)gNb(x,2,g));
    jumpT1 = mu*P0TangentJump(c4n,n4e,n4sDb,n4sNb,grad4e1,@(x)uDb(x,1,u4Db));
    jumpT2 = mu*P0TangentJump(c4n,n4e,n4sDb,n4sNb,grad4e2,@(x)uDb(x,2,u4Db));
    eta4s = length4s .* (jumpT1.^2+jumpT2.^2);%+jumpN1.^2+jumpN2.^2);

    %% Volume term
    f = f(mid4e);
    eta4e = area4e .* (f(:,1).^2+f(:,2).^2);

    %% Output 
    eta4e = sqrt(area4e.*eta4e + 1/2*sum(length4s(s4e).*eta4s(s4e),2));
end

function val = uDb(x,n,u4Db)
val = u4Db(x);
val = val(:,n);
end

function val = gNb(x,n,g)
val = g(x);
val = val(:,n);
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
