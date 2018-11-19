function afemP1PoissonLShapeExact
% afemP1PoissonLShapeExact.m
% Solve the Poisson equation on an L-Shape domain with linear P1 finite 
% elements adaptively.
%
% Given: function f. Seek for a  solution u such that
%   -div(grad(u)) = f   in Omega,
%              u  = 0   on Gamma_D,
%             u*n = g   on Gamma_N.
% with Dirichlet boundary Gamma_D and Neumann boundary Gamma_N. Compare the
% discrete solution with the exact solution u (in polar coordinates):
%               u = r^(2/3)*sin(2/3*phi)    
%% Initialization
    addpath(genpath(pwd));
    [c4n n4e n4sDb n4sNb] = loadGeometry('LshapeNb',1);
    minNrDoF = 5000;
    eta4nrDoF = sparse(1,1);

%% AFEM loop
    while( true )
        % SOLVE
        [x,nrDoF] = solveP1Poisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb); 
        %Exact error
        error4e = error4eP1L2(c4n,n4e,@uExact,x);
        error4nrDoF(nrDoF) = sqrt(sum(error4e));
		%Energy error
        en_error4e = error4eP1Energy(c4n,n4e,@gradExact,x);
        en_error4nrDoF(nrDoF) = sqrt(sum(en_error4e));
		% ESTIMATE
        [eta4s,n4s] = estimateP1EtaSides(@f,@g,@u4Db,x,c4n,n4e,n4sDb,n4sNb);
        eta4nrDoF(nrDoF) = sqrt(sum(eta4s));
        disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
            '; estimator = ',num2str(eta4nrDoF(nrDoF))]);
        if nrDoF >= minNrDoF, break, end;
        % MARK
        n4sMarked = markBulk(n4s,eta4s);
        % REFINE
        [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);
    end
    
%% Output error
    disp(['L2-Norm of the exact error: ',num2str(error4nrDoF(nrDoF))]); 
    disp(['L2-Norm of the energy error: ',num2str(en_error4nrDoF(nrDoF))]);
   
%% Plot solution and convergence graph.        
    figure;
    subplot(1,2,1);
    plotP1(c4n,n4e,uExact(c4n),'Exact Solution');
    subplot(1,2,2);
    plotP1(c4n,n4e,x,'P1 Solution');
    figure;
    plotP04e(c4n,n4e,error4e,'L2 error on elements');
    figure;
    plotP1(c4n,n4e,abs(uExact(c4n)-x),...
        'Difference of exact and approximated solution');
    figure;
    plotP04e(c4n,n4e,en_error4e,...
        {'Exact error of the gradients on elements',...
        '(grad u - grad u_l)'});
    nrDoF4lvl = find(eta4nrDoF);
    eta4lvl = eta4nrDoF(nrDoF4lvl);
    error4lvl=error4nrDoF(nrDoF4lvl);
    en_error4lvl = en_error4nrDoF(nrDoF4lvl);
    figure;
    plotConvergence(nrDoF4lvl,eta4lvl);
    hold all;
    plotConvergence(nrDoF4lvl,error4lvl);
    plotConvergence(nrDoF4lvl,en_error4lvl,'L2 energy error');

end

%% problem input data
function val = f(x)
   val = zeros(size(x,1),1);
end

function val = u4Db(x)
    [phi, r] = cart2pol(x(:,1),x(:,2));
    [index] = find(phi<0);
    phi(index) = phi(index) + 2*pi;
    val = r.^(2/3).*sin(2/3*phi);
end

function val = g(x)
    [phi, r] = cart2pol(x(:,1),x(:,2));
    [index] = find(phi<0);
    phi(index) = phi(index) + 2*pi;
    if ~isempty( find(phi<0)) && ~isempty(find(phi>2*pi))
        error('winkel')
    end
    for i = 1 : (size(x,1))
        if phi(i) >= 0 && phi(i) < 1/4*pi
            N = [1;0];   
        elseif phi(i) >= 1/4*pi && phi(i) < 3/4*pi
            N = [0;1];
        elseif phi(i) >= 3/4*pi && phi(i) < 5/4*pi
            N = [-1;0];
        elseif phi(i) >= 5/4*pi && phi(i) < 6/4*pi
            N = [0;-1];
        else
            error('rand')
        end
        N=[-sin(phi(i)), cos(phi(i)); cos(phi(i)), sin(phi(i))]*N;
        val(i,:) = 2/3*r(i)^(-1/3)*[  cos(2/3*phi(i)), ...
            sin(2/3*phi(i)) ]*N;               
    end
end

function val = uExact(x)
    [phi, r] = cart2pol(x(:,1),x(:,2));
    [index] = find(phi<0);
    phi(index) = phi(index) + 2*pi;
    val = r.^(2/3).*sin(2/3*phi);
end

function val = gradExact(x)
    [phi, r] = cart2pol(x(:,1),x(:,2));
    [index] = find(phi<0);
    phi(index) = phi(index)+2*pi;
    if ~isempty( find(phi<0)) && ~isempty(find(phi>2*pi))
        error('winkel')
    end
    val(:,1) = 2/3*r.^(-1/3).*cos(2/3*phi);
    val(:,2) = 2/3*r.^(-1/3).*sin(2/3*phi);
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
