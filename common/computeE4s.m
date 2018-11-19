function e4s = computeE4s(n4e)
%% computeE4s - Elements for sides.
%   computeE4s(n4e) returns a matrix each row of which corresponds to one side
%               of the decomposition. The side numbering is the same as in
%               n4s. Each row contains the numbers of the two elements that
%               the corresponding side is a part of. If it is a boundary
%               side the second entry is 0.
%               n4e is as specified in the documentation.
%
%   See also: computeN4s

    if isempty(n4e)
        e4s = zeros(0,2);
        return;
    end

    %% Compute e4s.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    [b,ind,back] = unique(sort(allSides,2),'rows','first');
    n4sInd = sort(ind); % by the way: n4s = allSides(n4sInd,:)

    nrElems = size(n4e,1);
    elemNumbers = [1:nrElems 1:nrElems 1:nrElems];
    e4s(:,1) = elemNumbers(n4sInd);
    allElem4s(ind) = accumarray(back,elemNumbers);
    e4s(:,2) = allElem4s(n4sInd)'-e4s(:,1);
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
