function [c4n n4e n4sDb n4sNb] = loadGeometry(name, OPTRefinementLevel)
%% loadGeometry - load data for a mesh.
% [c4n n4e n4sDb n4sNb] = loadGeometry('name') loads the data structures for
%   the mesh named 'name'. Optionally, the second parameter will
%   cause the mesh to be refined a given number of times using the uniform
%   red strategy.
% Example:
% [c4n n4e n4sDb n4sNb] = loadGeometry('LShape',2) loads the
%   mesh called 'LShape' and refines it two times.

    %% Load the geometry data.
    c4n = load([name,'_c4n.dat']);
    n4e = load([name,'_n4e.dat']);
    n4sDb = load([name,'_n4sDb.dat']);
    n4sNb = load([name,'_n4sNb.dat']);
    
    %% Initial refinement.
    if nargin < 2
        OPTRefinementLevel = 0;
    end
    for i=1:OPTRefinementLevel
        [c4n,n4e,n4sDb,n4sNb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
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
