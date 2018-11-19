function afemRT0Poisson
% afemRT0Poisson.m
% Solve the Poisson equation with RT0 P0 finite elements adaptively on
% a given geometry
%
% Seek for a  solution u such that
%   -div(grad(u)) = f   in Omega,
%              u  = 0   on Gamma_D,
%             u*n = g   on Gamma_N.
% with Dirichlet boundary Gamma_D and Neumann boundary Gamma_N.

    %% add paths
    addpath(genpath(pwd));

    close all;
    %% load the geometry
    %geom = 'Square';
    geom = 'Lshape';
    %geom = 'Slit';
    [c4n, n4e, n4sDb, n4sNb] = loadGeometry(geom,1);

    %% set the minimal number of nodes
    minNrDoF = 1000;

    %% initialisation
    nrDoF4lvl = [];
    eta4lvl = [];

    %% AFEM loop
    tic
    while( true )
        % SOLVE
        [p,u,nrDoF] = solveRT0Poisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
        nrDoF4lvl(end+1) = nrDoF;
        % ESTIMATE
        [eta4s,n4s] = estimateRT0EtaSides(@f,@g,@u4Db,p,u,c4n,n4e,n4sDb,n4sNb);
        eta4lvl(end+1) = sqrt(sum(eta4s));
        disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
              '; estimator = ',num2str(eta4lvl(end))]);
        % leave the loop if minNrDoF is reached
        if nrDoF >= minNrDoF, break, end;
        % MARK
        n4sMarked = markBulk(n4s,eta4s);
        % REFINE
        [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);
    end

    toc

    %% plot
    figure;
    plotTriangulation(c4n,n4e);
    figure;
    plotP04e(c4n,n4e,u,...
       {'RT0 Solution - u'...
       [num2str(length(p) + length(u)) ' degrees of freedom']});
    figure;
    plotRT04e(c4n,n4e,p,{'RT0 Solution - p'...
        [num2str(length(p) + length(u)) ' degrees of freedom']});
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
    val = zeros(size(x,1),1);
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
