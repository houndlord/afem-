function length4s = computeLength4s(c4n,n4s)
%% computeLength4s - Lengths for sides.
%   computeLength4s(c4n, n4s) computes the length of each side of the
%                         decomposition. c4n and n4s are as specified in the
%                         documentation.
%
%   See also: computeN4s, computeArea4e
    
    if isempty(n4s)
        length4s = zeros(0,1);
        return;
    end

    %% Compute length4s in a vectorised manner.
    length4s = sqrt( sum( (c4n(n4s(:,2),:) - c4n(n4s(:,1),:)).^2, 2) );
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
