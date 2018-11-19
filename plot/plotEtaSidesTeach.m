function plotEtaSidesTeach(figure, c4n, n4s, eta4s)
%% Draw the estimated Error on sides.
%   plotEtaSides(figure, c4n, n4s, eta4s) draws the estimated Error defined
%                       by the grid (c4n, n4s) and the error for each side
%                       into the given figure.

    angle = [-28, 56];  % set the point of view
    
    clf;
    [errX,errY,errZ,errC] = getErrorXYZC(c4n, n4s, eta4s);
    patch(errX,errY,errZ,errC);

    set(figure,'Name','AFEM-Error','Position',[670 500 600 400]);
    title('Estimated error on sides');
    view(angle);
    drawnow;
end

%% Function to get the coordinates for the patch-function
function [valX, valY, valZ,valC] = getErrorXYZC(c4n, n4s, eta4s)
    % coordinates for all nodes of sides
    X1 = c4n(n4s(:,1), 1);
    X2 = c4n(n4s(:,2), 1);
    Y1 = c4n(n4s(:,1), 2);
    Y2 = c4n(n4s(:,2), 2);   

    % sides for the error in patch-style
    valX = [X1';X1';X2';X2'];
    valY = [Y1';Y1';Y2';Y2'];

    % each side should have the same height - the error
    valZ = [zeros(size(eta4s'));eta4s';eta4s';zeros(size(eta4s'))];
    valC = valZ / max(eta4s);
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
