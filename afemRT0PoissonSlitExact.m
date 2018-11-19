function afemRT0PoissonSlitExact
% afemRT0Poisson.m
% Solve the Poisson equation with RT0 P0 finite elements adaptively on
% a given geometry
%
% Seek for a  solution u such that
%   -div(grad(u)) = f   in Omega,
%              u  = 0   on Gamma_D,
%             u*n = g   on Gamma_N.
% with Dirichlet boundary Gamma_D and Neumann boundary Gamma_N. Compare the
% discrete solution with the exact solution u (in polar coordinates):
%               u = r^(1/4)*sin(1/4*phi)

%% Initialization
    addpath(genpath(pwd));
    close all;
    [c4n, n4e, n4sDb, n4sNb] = loadGeometry('SlitNb',1);
    minNrDoF = 1000;
    nrDoF4lvl = [];
    eta4lvl = [];
    error4lvl = [];
    energy4lvl = [];

%% AFEM loop
    tic
    while( true )
        % SOLVE
        [p,u,nrDoF] = solveRT0Poisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
        nrDoF4lvl(end+1) = nrDoF;
        %Exact error
        error4e = error4eRT0L2(c4n, n4e, @uExact, u);
        error4lvl(end+1) = sqrt(sum(error4e));
        %Energy error
        energy4e = error4eRT0Energy(c4n, n4e, @gradExact, p);
        energy4lvl(end+1) = sqrt(sum(energy4e));
        % ESTIMATE
        [eta4s,n4s] = estimateRT0EtaSides(@f,@g,@u4Db,p,u,c4n,n4e,n4sDb,n4sNb);
        eta4lvl(end+1) = sqrt(sum(eta4s));
        disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
              '; estimator = ',num2str(eta4lvl(end))]);
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
    hold all;
    plotConvergence(nrDoF4lvl,error4lvl,'||u - u_l||_{L2}');
    plotConvergence(nrDoF4lvl,energy4lvl,'||\nablau - \nablau_l||_{L2}');
end

%% problem input data
function val = f(x)
    [phi, r] = cart2pol(x(:,1),x(:,2));
    phi(phi < -eps) = phi(phi < -eps)+2*pi;
    % laplace u = d^2u/dr^2 + 1/r*du/dr + 1/r^2 * d^2u/dphi^2
    val = zeros(size(x,1),1);
end


function val = u4Db(x)
    val = uExact(x);
end


function val = g(x)
    [phi,r] = cart2pol(x(:,1),x(:,2));
    phi(phi < -eps) = phi(phi < -eps)+2*pi;
    val = zeros(size(x,1),1);
    for i = 1 : (size(x,1))
        if (x(i,1)==1)
            N = [1;0];
        elseif (x(i,2)==1)
            N = [0;1];
        elseif (x(i,1)==-1)
            N = [-1;0];
        elseif (x(i,2)==-1)
            N = [0;-1];
        elseif (x(i,2)==0)
            N = [0;-1];
        else
            error('Normalen am Neumannrand');
        end
        val(i) = gradExact(x) * N;
    end
end


function val = uExact(x)
    [phi, r] = cart2pol(x(:,1),x(:,2));
    phi(phi < -eps) = phi(phi < -eps)+2*pi;
    val = r.^(1/4).*sin(1/4*phi);
end


function val = gradExact(x)
    x1 = x(:,1);
    x2 = x(:,2);
    [phi, r] = cart2pol(x1,x2);
    phi(phi < -eps) = phi(phi < -eps)+2*pi;
    if (~isempty(phi(phi<-eps)) | ~isempty(phi(phi>2*pi)))
        error('umrechnung in polarkoordinaten')
    end
    val(:,1) = 1/4*r.^(-3/4).*cos(1/4*phi);
    val(:,2) = 1/4*r.^(-3/4).*sin(1/4*phi);
    for i = 1 : (size(x,1))
        if r(i)==0
            val(i,2) = Inf;
        end
        val(i,:)=val(i,:)*[-sin(phi(i)), cos(phi(i));...
            cos(phi(i)), sin(phi(i))];
    end
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
