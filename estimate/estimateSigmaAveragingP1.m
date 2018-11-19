function eta4e = estimateSigmaAveragingP1(c4n,n4e,sigma4e)
%% estimateSigmaAveragingP1 - L2 difference of sigma and A(sigma)
% Given an elementwise constant function sigma on a triangulation 
% ompute the L2 norm of A(sigma)-sigma.
% Input: c4n,n4e: triangulation
%        sigma4e: values of sigma for each element
%
% Output: eta4e: ||A(sigma)-sigma||_L2(T) for each element T

    %% Initialisation
    nrElems = size(n4e,1);
    area4e  = computeArea4e(c4n,n4e);
    eta4e   = zeros(nrElems,1);
    
    %% Compute the P1 average of sigma
    A = P0AveragingP1(c4n,n4e,sigma4e);

    %% Eta for each element
    for elem = 1 : nrElems
        s=0.5*[1 1 0;0 1 1;1 0 1]*A(n4e(elem,:),:)-[1;1;1]*sigma4e(elem,:);
        eta4e(elem) = sum(sum(s.^2))*area4e(elem)/3;
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
