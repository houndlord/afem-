function s4e = computeS4e(n4e)
%% computeS4e - Sides for elements.
%   computeS4e(n4e) returns a matrix each row of which corresponds to one
%               element of the decomposition. Each row contains the numbers
%               of the three sides belonging to an element. The side
%               numbering is the same as in n4s.
%
%   See also: computeN4s

    if isempty(n4e)
        s4e = [];
        return;
    end

    %% Compute s4e.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    [b,ind,back] = unique(sort(allSides,2),'rows','first');
    [n4sInd, sortInd] = sort(ind); % by the way: n4s = allSides(n4sInd,:)
    sideNr(sortInd) = 1:length(ind); % sideNr(back): numbers for allSides
    s4e = reshape(sideNr(back),size(n4e));
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
