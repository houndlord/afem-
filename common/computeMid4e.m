function mid4e = computeMid4e(c4n, n4e)
%% computeMid4e - midpoints for elements.
%   computeMid4e(c4n, n4e) computes the midpoint for each element of the
%                     decomposition. c4n and n4e are as specified in the
%                     documentation.
%
%   See also: computeArea4e, computeMid4s

    if isempty(n4e)
        mid4e = zeros(0,2);
        return;
    end

    %% Compute mp4e.
    mid4e = ( c4n(n4e(:,1),:) + c4n(n4e(:,2),:) + c4n(n4e(:,3),:) ) / 3;
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
