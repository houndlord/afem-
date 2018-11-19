function plotP04e(c4n, n4e, x, OPTtitle)
%% Draw the P0-function on elements.
%   plotP04e(c4n, n4e, x, OPTtitle) draws the P0-function defined by
%                                   the grid (c4n, n4e) and the coefficient
%                                   vector x into the given figure. The
%                                   input argument OPTtitle is optional, it
%                                   sets the title of the figure. The
%                                   default value is empty.

    angle = [-37.5, 30];

    [errX,errY,errZ,errC] = getErrorXYZC(c4n, n4e, x);
    if(size(n4e,1)>2000)
        patch(errX,errY,errZ,errC,'EdgeColor','none');
    else
        patch(errX,errY,errZ,errC);
    end


    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end

    view(angle);
    grid on;
    drawnow;
end

%% Function to get the coordinates for the patch-function
function [valX, valY, valZ,valC] = getErrorXYZC(c4n, n4e, x)
    % coordinates for all nodes of sides
    X1 = c4n(n4e(:,1), 1);
    X2 = c4n(n4e(:,2), 1);
    X3 = c4n(n4e(:,3), 1);
    Y1 = c4n(n4e(:,1), 2);
    Y2 = c4n(n4e(:,2), 2);
    Y3 = c4n(n4e(:,3), 2);

    % triangles for the function in patch-style
    valX = [X1';X2';X3';X3'];
    valY = [Y1';Y2';Y3';Y3'];

    % each element should have constant height - the function value
    valZ = [x';x';x';x'];
    valC = valZ / max(x);

    % create walls.
    nrElems = size(n4e,1);
    n4s= [n4e(:,1:2);n4e(:,2:3); n4e(:,3),n4e(:,1)];
    X  = reshape(c4n(n4s',1),2,3*nrElems);
    Y  = reshape(c4n(n4s',2),2,3*nrElems);
    X  = [X; X(2,:); X(1,:)];
    Y  = [Y; Y(2,:); Y(1,:)];
    Z  = [1; 1; 0; 0]*[x', x', x'];

    valX = [valX, X];
    valY = [valY, Y];
    valZ = [valZ, Z];
    valC = [valC, Z/max(x)];
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
