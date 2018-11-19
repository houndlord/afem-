function afemCRPoisson
%% afemCRPoisson - Solve Poisson model problem with CR elements.
% Solve the Poisson equation with linear CR finite elements adaptively on
% the domain Omega.
%
% Seek for a  solution u such that
%   -div(grad(u)) = f   in Omega,
%              u  = 0   on Gamma_D,
%             u*n = g   on Gamma_N.
% with Dirichlet boundary Gamma_D and Neumann boundary Gamma_N.

    %% Initialization
    addpath(genpath(pwd));
    [c4n, n4e, n4sDb, n4sNb] = loadGeometry('Lshape',1);
    minNrDoF = 1000;
    eta4nrDoF = sparse(1,1);

    %% AFEM loop
    while( true )
        % SOLVE
        [x,nrDoF] = solveCRPoisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
        % ESTIMATE
        [eta4s,n4s] = estimateCREtaSides(@f,@g,@u4Db,x,c4n,n4e,n4sDb,n4sNb);
        eta4nrDoF(nrDoF) = sqrt(sum(eta4s));
        disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
            '; estimator = ',num2str(eta4nrDoF(nrDoF))]);
        if nrDoF >= minNrDoF, break, end;
        % MARK
        n4sMarked = markBulk(n4s,eta4s);
        % REFINE
        [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);
    end

    %% Plot mesh, solution and convergence graph.
    figure;
    plotTriangulation(c4n,n4e);
    figure;
    plotCR(c4n,n4e,x,{'CR Solution'; [num2str(nrDoF) ' degrees of freedom']});
    nrDoF4lvl = find(eta4nrDoF);
    eta4lvl = eta4nrDoF(nrDoF4lvl);
    figure;
    plotConvergence(nrDoF4lvl,eta4lvl,'\eta_l');
end

%% problem input data
function val = f(x)
    val = ones(size(x,1),1);
end

function val = u4Db(x)
    val = zeros(size(x,1),1);
end

function val = g(x)
    val = ones(size(x,1),1);
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
