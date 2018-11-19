function afemP1PoissonShort

addpath(genpath(pwd));
[c4n, n4e, n4sDb, n4sNb] = loadGeometry('Square',1);
eta4nrDoF = sparse(1,1);

for l = 1 : 6
    % SOLVE
    [x,nrDoF] = solveP1Poisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
    % ESTIMATE
    [eta4s,n4s] = estimateP1EtaSides(@f,@g,@u4Db,x,c4n,n4e,n4sDb,n4sNb);
    eta4nrDoF(nrDoF) = norm(eta4s);
    % display
    disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
        '; estimator = ',num2str(eta4nrDoF(nrDoF))]);
    figure;
    plotP1(c4n,n4e,x,{'P1-Solution' [num2str(nrDoF) ' degrees of freedom']});pause(0.1)
    % MARK
    n4sMarked = markBulk(n4s,eta4s);
    % REFINE
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);

end

figure;
plotConvergence(find(eta4nrDoF),nonzeros(eta4nrDoF)','\eta_l');
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
