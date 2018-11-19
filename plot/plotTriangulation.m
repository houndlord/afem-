function plotTriangulation ( c4n, n4e )
%% Draw a triangular grid into a new figure.
%   plotTriangulation(c4n, n4e) draws the grid defined by c4n and n4e into
%                               a figure.

    % Set titles of plot and window.
    title({'Mesh plot'; [num2str(size(c4n,1)), ' nodes']});
    % This can be done with triplot but patch is _much_ faster.
    % Get the coordinates for each node of each triangle.
    X1 = c4n(n4e(:,1),1);
    Y1 = c4n(n4e(:,1),2);
    X2 = c4n(n4e(:,2),1);
    Y2 = c4n(n4e(:,2),2);
    X3 = c4n(n4e(:,3),1);
    Y3 = c4n(n4e(:,3),2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];

    % Set the colour each triangle is filled with.
    C = 'white';

    % Draw everything, make sides blue (looks more like triplot).
    patch(X,Y,C,'EdgeColor','blue');
    axis equal tight;
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
