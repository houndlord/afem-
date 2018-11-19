function s4n = computeS4n(n4e, n4s)
%% computeS4n - Sides for nodes.
%   computeS4n(n4e) returns a symmetric sparse matrix in which the entry (j,k)
%               contains the number of the side with the end nodes j and k
%               or zero if no such side exists.
%               The side numbering is the same as in n4s. n4e is as
%               specified in the documentation.
%
%   See also: computeN4s, computeS4e

  if isempty(n4e)
      s4n = [];
      return;
  end

  %% Optionally compute n4s.
  if nargin < 2
      n4s = computeN4s(n4e);
  end

  %% Compute s4n.
  S = size(n4s, 1);
  N = max(n4e(:));
  s4n = sparse(n4s(:, 1), n4s(:, 2), 1:S, N, N);
  % Up to here, s4n is not yet symmetric as each side has only
  % been considered once. The following makes sure that s4n is
  % symmetric.
  s4n = s4n + s4n';

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
