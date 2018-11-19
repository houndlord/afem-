function plotCR(c4n, n4e, x, OPTtitle)
%% Draw a Crouzeix-Raviart-function.
%   plotCR(c4n, n4e, x, OPTtitle) draws the CR-function defined by the grid
%                                 (c4n, n4e) and the basis coefficients (x)
%                                 The input argument OPTtitle is optional,
%                                 it sets the title of the figure. The
%                                 default value is empty.

    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end

    %% Get coordinates for nodes.
    X1 = c4n(n4e(:,3),1);
    Y1 = c4n(n4e(:,3),2);
    X2 = c4n(n4e(:,1),1);
    Y2 = c4n(n4e(:,1),2);
    X3 = c4n(n4e(:,2),1);
    Y3 = c4n(n4e(:,2),2);

    %% Translate values for degrees of freedom into values for nodes.
    s4e = computeS4e(n4e);
    W = x(s4e)';    % Get x for each side of each element.
    Z = (ones(3)-2*eye(3))*W;

    %% Assemble parameters for the patch function.
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    % The colour of a triangle is determined by its midpoint.
    C = sum(Z,1)/3;

    %% Put everything together.
    % For large numbers of elements, omit the black boundary around each
    % triangle.
    if( size(n4e,1) > 2000 )
        patch(X,Y,Z,C,'EdgeColor','none');
    else
        patch(X,Y,Z,C);
    end
    view(-37.5,30);
    grid on;
    drawnow;
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
