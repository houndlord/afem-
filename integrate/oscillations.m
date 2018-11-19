function [osc4e ,mean4e] = oscillations(c4n,n4e,f,degree)
%% oscillcations - oscillations of a function on a mesh
%  Input:    c4n,n4e - mesh
%            f       - R^2->R; input: points; output: values
%            degree  - accuracy of integration
%  Output:   osc4e   - vector of oscillations squared for each element
%	     mean4e  - vector of integral means of f for each element
%% Initialisation
    meanDeg = degree;
    oscDeg  = 2*degree;
    area4e  = computeArea4e(c4n,n4e);
%% Compute the integral mean of f on each element.
    mean4e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)(f(Gpts4p)),meanDeg)...
             ./area4e;
%% Compute oscillations: locally on each T - osc4e = ||f-mean(f)||_L2(T) / |T|
    osc4e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)((f(Gpts4p)-mean4e).^2),...
    	     oscDeg);
    osc4e = osc4e .* area4e;
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
