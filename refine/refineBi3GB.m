function [c4nNew,n4eNew,n4sDbNew,n4sNbNew] = ...
			            refineBi3GB(c4n,n4e,n4sDb,n4sNb,n4sMarked)
%% refineBi3GB - refine using the Bisec3-Green-Blue-strategy
%   refineBi3GB(c4n, n4e, n4sDb, n4sNb, n4sMarked) Refines a given mesh using
%       the Bisec3-Green-Blue refinement. For details on data structures
%       and refinement strategies see the documentation. Input is a mesh
%       defined by c4n, n4e, n4sDb, n4sNb and marked sides given by
%       n4sMarked. Output is a refined mesh defined by c4nNew, n4eNew,
%       n4sDbNew and n4sNbNew.

    nrNodes = size(c4n,1);
    nrElems = size(n4e,1);

    %% Closure
    n4sRefine = closure(n4e,n4sMarked);

    %% Compute newNodes4n to find new nodes faster.
    newNodes4s = sparse(n4sRefine(:,1),n4sRefine(:,2),...
                       (1:size(n4sRefine,1))'+ nrNodes, nrNodes,nrNodes);
    newNodes4s = newNodes4s + newNodes4s';

    %% Compute coordinates of new nodes.
    mid4sRefine = (c4n(n4sRefine(:,1),:)+c4n(n4sRefine(:,2),:))/2;
    c4nNew = [c4n;mid4sRefine];

    %% bisec3 refinement
    % Count elements in new triangulation. For each refined inner side, two
    % new elements will be created.
    nrNewElems = nrElems+2*size(n4sRefine,1)-size(n4sDb,1)-size(n4sNb,1);
    n4eNew = zeros(nrNewElems,3);
    ind = 0; % index to keep track of the current element number in n4eNew
    for curElem = 1 : nrElems
        curNodes = n4e(curElem,:);
        curNewNodes = [newNodes4s(curNodes(1),curNodes(2));
                       newNodes4s(curNodes(2),curNodes(3));
                       newNodes4s(curNodes(3),curNodes(1));
                      ];
	nrNewNodes4curElem = nnz(curNewNodes);
        if nrNewNodes4curElem == 0 % no refinement
            n4eNew(ind+1,:) = curNodes;
            ind = ind+1;
        elseif nrNewNodes4curElem == 1 % green refinement
            n4eNew(ind+1:ind+2,:) = ...
                [ curNodes(3)    curNodes(1) curNewNodes(1);
                  curNodes(2)    curNodes(3) curNewNodes(1);
                ];
            ind = ind+2;
        elseif nrNewNodes4curElem == 2
            if curNewNodes(2) > 0 % blue right
                n4eNew(ind+1:ind+3,:) = ...
                    [ curNodes(3)    curNodes(1)    curNewNodes(1);
                      curNewNodes(1) curNodes(2)    curNewNodes(2);
                      curNodes(3)    curNewNodes(1) curNewNodes(2);
                    ];
            else % blue left
                n4eNew(ind+1:ind+3,:) = ...
                    [ curNodes(1)    curNewNodes(1) curNewNodes(3);
                      curNewNodes(1) curNodes(3)    curNewNodes(3);
                      curNodes(2)    curNodes(3)    curNewNodes(1);
                    ];
            end
            ind = ind+3;
        elseif nrNewNodes4curElem == 3 % bisec3 refinement
            n4eNew(ind+1:ind+4,:) = ...
                [ curNewNodes(1) curNodes(2)    curNewNodes(2);
                  curNodes(3)    curNewNodes(1) curNewNodes(2);
                  curNodes(1)    curNewNodes(1) curNewNodes(3);
                  curNewNodes(1) curNodes(3)    curNewNodes(3);
                ];
            ind = ind+4;
        end
    end

    %% refinement of  Dirichlet boundary
    nrNewDbSides= size(intersect(sort(n4sDb,2),sort(n4sRefine,2),'rows'),1);
    n4sDbNew = zeros(nrNewDbSides,2);
    ind = 0;
    for curSide = 1 : size(n4sDb,1)
        curNodes = n4sDb(curSide,:);
        curNewNodes = newNodes4s(curNodes(1),curNodes(2));
        if curNewNodes == 0
            n4sDbNew(ind + 1,:) = curNodes;
            ind = ind + 1;
        else
            n4sDbNew(ind + 1:ind + 2,:) = ...
                [ curNodes(1)    curNewNodes;
                  curNewNodes    curNodes(2);
                ];
            ind = ind + 2;
        end
    end

    %% refinement of  Neumann boundary
    nrNewNbSides= size(intersect(sort(n4sNb,2),sort(n4sRefine,2),'rows'),1);
    n4sNbNew = zeros(nrNewNbSides,2);
    ind = 0;
    for curSide = 1 : size(n4sNb,1)
        curNodes = n4sNb(curSide,:);
        curNewNodes = newNodes4s(curNodes(1),curNodes(2));
        if curNewNodes == 0
            n4sNbNew(ind + 1,:) = curNodes;
            ind = ind + 1;
        else
            n4sNbNew(ind + 1:ind + 2,:) = ...
                [ curNodes(1)    curNewNodes;
                  curNewNodes    curNodes(2);
                ];
            ind = ind + 2;
        end
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
