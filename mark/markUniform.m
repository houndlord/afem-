function n4sMarked = markUniform(n4p)
%% markUniform - Mark given parts uniformly.
%   n4sMarked = markUniform(n4p) marks every given part for refinement,
%       regardless of estimated errors. n4p contains the nodes for the
%       parts: 2 entries per row for sides, 3 entries per row for elements.
%       The output is a list of marked sides given by their end nodes.

    nrParts = size(n4p,1);
    dimParts = size(n4p,2);

    %% Uniform criterion: mark everything
    I = 1:nrParts;

    %% Mark sides
    if dimParts == 2 % Sides were given. Mark 'em all.
        n4sMarked = n4p(I,:);
    elseif dimParts == 3 % Elements were given. Mark all sides of these.
        allSidesMarked = [n4p(I,[1 2]);n4p(I,[2 3]);n4p(I,[3 1])];
        % Eliminate duplicates.
        [b, ind]   = unique(sort(allSidesMarked,2), 'rows');
        n4sMarked = allSidesMarked(sort(ind),:);
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
