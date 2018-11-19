function val = integrate(c4n, n4p, integrand, degree, OPTsize4parts)
% Integrates function handle integrand over given domain 1D or 2D.
% 'integrand' must have the form integrand(n4p,Gpts4p,Gpts4ref) and the
% return value must have the dimension [ nrParts n m ]
%
% integrate.m
% input:  c4n           - coordinates for nodes
%         n4p           - nodes of the parts of a partition of a 1D or 2D domain
%         integrand     - function to integrate
%         degree        - integration is exact for polynomials up to degree
%         OPTsize4parts - length4s (1D) or area4e (2D) (Optional parameter)
% output: val           - integral value of dimension [ nrParts n m ]

% [val of dimension (nrParts x n x m)] = integrand(n4p,Gpts4p,Gpts4ref)
% input: n4p            - nodes of the parts of a partition of a 1D or 2D domain
%        Gpts4p         - coordinates of the transformed GaussPoints on each part
%        Gpts4ref       - coordinates of the GaussPoints of the reference intervall
%                         or the reference triangle

    %% Number of Parts
    nrParts = size(n4p,1);

    %% Number of Gauss Points 
    nrGpts = ceil((degree+1)/2);

    %% 1D or 2D 
    if size(n4p,2) == 2 
        % 1D
        [Gpts4ref,Gwts4ref] = getGaussPoints(nrGpts);
        if nargin == 5
            weightFactor = OPTsize4parts;
        else
            weightFactor = computeLength4s(c4n,n4p);
        end

    elseif size(n4p,2) == 3 
        % 2D triangle
        [Gpts4ref,Gwts4ref] = getConProdGaussPoints(nrGpts);
        if nargin == 5 
            weightFactor = 2*OPTsize4parts;
        else
            weightFactor = 2*computeArea4e(c4n,n4p);
        end
        nrGpts = nrGpts^2;

    else
        error('Could not integrate because no Domain is specified.');
    end
    % Integrate 
    for curGpt = 1:nrGpts
        %transformation
        curGpt4p = ref2arbitrary(c4n,n4p,Gpts4ref(curGpt,:));
        % evaluate function handle eval must be [nrParts n m]
        integrandVal = integrand(n4p,curGpt4p,Gpts4ref(curGpt,:));
        if(curGpt == 1)
            dim = size(integrandVal);
            if numel(dim) < 3; dim(3) = 1; end
            val = zeros(nrParts,dim(2),dim(3));
        end
        % weighted sum
        curGwt4p = weightFactor*Gwts4ref(curGpt);
        curGwt4p = reshape( curGwt4p*ones(1,dim(2)*dim(3)),[nrParts,dim(2),dim(3)]);
        val = val + curGwt4p.*integrandVal;
    end
end


%% ref2arbitrary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = ref2arbitrary(c4n,n4p,curGpt)
    %1D or 2D
    if size(n4p,2) == 2
        %1D ref = conv(0,1)
        c1 = c4n(n4p(:,1),:);
        c2 = c4n(n4p(:,2),:);

        val = c1 + curGpt(1)*(c2-c1);   

    elseif size(n4p,2) == 3 
        % 2D ref = conv{ (0,0),(1,0),(0,1)}
        c1 = c4n(n4p(:,1),:);
        c2 = c4n(n4p(:,2),:);
        c3 = c4n(n4p(:,3),:);

        val = c1 + curGpt(1)*(c2-c1) + ...
                   curGpt(2)*(c3-c1);
    end
end

function [x,w] = getGaussPoints(n)
% Gives n Gauss points for the unite vector.
% Used for Gauss-Legendere integration.
%
% input:  n - number of points
% output: x - Gauss points
%         w - Gauss weights

    % Find n Gauss_Legendre Points for Intervall [-1,1]
    gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
    [V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
    x = diag(D);
    w = 2*V(1,:).^2;

    % linear map to Intervall [0,1]
    x = .5 * x + .5;
    w = .5 * w';
end

function [x,w] = getConProdGaussPoints(n)
% Gauss Points for the Reference Triangle
%    T = {x+y| 0<=x<=1, 0<=y<=1, x+y<=1}
% Integrates polynomials up to degree 2n-1 exact
% using the Stroud Conical Product rule with n^2
% quadrature Points.
%
% input:  n - number of points
% output: x - Gauss points
%         w - Gauss weights

    % Find n Gauss_Legendre Points for Intervall [-1,1]
    gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
    [V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
    r = diag(D);
    a = 2*V(1,:).^2;% norm-factor -> int(1,-1,1)=2

    % Find n Gauss_Jacobi Points for Intervall [-1,1]
    delta = -1./(4*(1 : n).^2-ones(1,n));
    gamma = sqrt((2 : n).*(1 : n-1)) ./ (2*(2 : n)-ones(1,n-1));
    [V,D] = eig( diag(delta)+diag(gamma,1)+diag(gamma,-1) );
    s = diag(D);
    b = 2*V(1,:).^2; % norm-factor -> int((1-x),-1,1)=2

    % linear map to Intervall [0,1]
    % w(x) = 1 changes norm-factor from 2 to 1
    r = .5 * r + .5;
    a = .5 * a';
    % w(x) = (1-x) changes norm-factor from 2 to 1/2
    s = .5 * s + .5;
    b = .25 * b';

    % conical product [ s_j , r_i(1-s_j) ]  a_i*b_j
    s = repmat(s',n,1); s = s(:);
    r = repmat(r,n,1);
    x = [ s , r.*(ones(n^2,1)-s) ];
    w = a*b';
    w = w(:);
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
