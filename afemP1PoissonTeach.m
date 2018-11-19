function afemP1PoissonTeach
% afemP1PoissonTeach.m
% Solve the Poisson equation with linear P1 finite elements adaptively on the.
%
% Seek for a  solution u such that
%   -div(grad(u)) = f   in Omega,
%              u  = 0   on Gamma_D,
%             u*n = g   on Gamma_N.
% with Dirichlet boundary Gamma_D and Neumann boundary Gamma_N.
%
% In addition the steps are visulized by 5 plots per level: the
% triangulation, the error estimates, the marked sides, the marked sides
% after the closure algorithm and the colored triangles (RGB)

    %% add paths
    addpath(genpath(pwd));

    %% load the geometry
    % geom = 'Square';
    geom = 'Lshape';
    % geom = 'Slit';
    [c4n, n4e, n4sDb, n4sNb] = loadGeometry(geom,1);

    %% set the maximal number of nodes
    minNrDoFs = 1000;

    %% initialisation
    nrDoF4lvl = [];
    eta4lvl = [];
    
    %% AFEM loop
    % Solve and estimate at least once. Decide whether or not to continue
    % the AFEM loop afterwards.
    tic   

    % solve
    [x,nrDoF] = solveP1Poisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
    nrDoF4lvl(end+1) = nrDoF;
    % plot first triangulation                 
    aFemLoopFigure = figure;
    set(aFemLoopFigure,'Name','AFEM-Teach','Position',[50 500 600 400]);
    angle = [-28, 56];
    [triX,triY,triZ] = getTriangulationXYZ(c4n, n4e);       
    patch(triX,triY,[1 1 1]); %white triangles 
    myaxis = axis;
    
    % estimate
    [eta4s,n4s] = estimateP1EtaSides(@f,@g,@u4Db,x,c4n,n4e,n4sDb,n4sNb);
    eta4lvl(end+1) = norm(eta4s);
    disp(['nodes/dofs: ',int2str(size(c4n,1)),'/',num2str(nrDoF),...
          '; estimator = ',num2str(eta4lvl(end))]);
    % plot first error estimates    
    aFemErrorFigure = figure;
    plotEtaSidesTeach(aFemErrorFigure, c4n, n4s, eta4s);
      
    % While we have not reached the desired number of degrees of freedom
    % yet, execute the AFEM loop.
    while( nrDoF < minNrDoFs )

        % mark
        n4sMarked = markBulk(n4s,eta4s);
        % plot marked sides
        waitForUser();                
        figure(aFemLoopFigure);
        set(0,'CurrentFigure',aFemLoopFigure);
        clf;                                                                  
        [triX,triY,triZ] = getTriangulationXYZ(c4n, n4e);                    
        patch(triX,triY,triZ,[1 1 1]);   
        [markX,markY,markC] = getMarkXYC(c4n, n4sMarked);   
        axis(myaxis);
        patch(markX,markY,markC,'EdgeColor','none');
        
        % plot reference sides
        waitForUser();
        figure(aFemLoopFigure);
        set(0,'CurrentFigure',aFemLoopFigure);
        [refX,refY,refC] = getRefXYC(c4n, n4e);  
        axis(myaxis);
        patch(refX,refY,refC);        
        
        % plot marked sides after closure --> new edges with other color
        n4sRefine = closure(n4e,n4sMarked);              % just for the plot!                      
        waitForUser();     
        figure(aFemLoopFigure);
        set(0,'CurrentFigure',aFemLoopFigure);
        [markX,markY,markC] = getMarkXYC(c4n, n4sRefine);     
        axis(myaxis);    
        markC(1,:,:) = markC(1,:,[2 3 1]);               % changes color to blue
        [n4sClosure, indexN4sRefine, indexN4sMarked] = intersect(n4sRefine, [n4sMarked; n4sMarked(:,[2 1])], 'rows');
        markC(1,indexN4sRefine,:) = ones(length(indexN4sRefine),1)*[1 0 0];    % already marked sides red again
        patch(markX,markY,markC,'EdgeColor','none');
        
        % plot how afem will refine
        waitForUser();      
        figure(aFemLoopFigure);
        set(0,'CurrentFigure',aFemLoopFigure);
        clf;        
        [colX,colY,colC] = getColoredXYC(c4n, n4e, n4sRefine);             
        patch(colX,colY,colC);         
        %refine
        [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);  
        
        % plot new triangulation
        waitForUser();  
        figure(aFemLoopFigure);
        set(0,'CurrentFigure',aFemLoopFigure);
        clf;
        [triX,triY,triZ] = getTriangulationXYZ(c4n, n4e);                            
        patch(triX,triY,triZ,[1,1,1]); 
                
        % solve
        [x,nrDoF] = solveP1Poisson(@f,@g,@u4Db,c4n,n4e,n4sDb,n4sNb);
        nrDoF4lvl(end+1) = nrDoF;            
        % estimate
        [eta4s,n4s] = estimateP1EtaSides(@f,@g,@u4Db,x,c4n,n4e,n4sDb,n4sNb);
        eta4lvl(end+1) = norm(eta4s);
        disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
              '; estimator = ',num2str(eta4lvl(end))]);
        % plot estimated error 
        figure(aFemErrorFigure);
        set(0,'CurrentFigure',aFemErrorFigure);
        plotEtaSidesTeach(aFemErrorFigure, c4n, n4s, eta4s);
    end

    toc    
    
    figure;   
    plotP1(c4n,n4e,x,{'P1-Solution' [num2str(length(x)) ' nodes']});
    figure;
    plotConvergence(nrDoF4lvl,eta4lvl,'Eta');
end

%% problem input data
function val = f(x)
    val = ones(size(x,1),1);
end

function val = u4Db(x)
    val = zeros(size(x,1),1);
end

function val = g(x)
    val = zeros(size(x,1),1);
end

%% functions for teach-plot 
% function which returns the input parameters of the patch function (no Color) to
% draw a triangulation - also used as basis for the marking
function [valX, valY, valZ] = getTriangulationXYZ(c4n, n4e)
    % coordinates for triangles 
    X1 = c4n(n4e(:,1), 1);
    X2 = c4n(n4e(:,2), 1);
    X3 = c4n(n4e(:,3), 1);
    Y1 = c4n(n4e(:,1), 2);
    Y2 = c4n(n4e(:,2), 2); 
    Y3 = c4n(n4e(:,3), 2);     
    
    % combined coordinates in patch-style
    valX = [X1';X2';X3'];
    valY = [Y1';Y2';Y3'];
    % no heigth
    valZ = zeros(size(valX));
end

% function which returns the input parameters of the patch function to
% mark the reference sides with a second parallel line.
% reference sides are the side between the first two nodes in n4e
function [valX, valY, valC]  = getRefXYC(c4n, n4e)
    % coordinates of reference sides
    X1 = c4n(n4e(:,1), 1);
    X2 = c4n(n4e(:,2), 1);
    Y1 = c4n(n4e(:,1), 2);
    Y2 = c4n(n4e(:,2), 2);     
    
    % combined coordinates in patch-style
    valX = [X1';X2'];
    valY = [Y1';Y2'];
        
    for i=1 : size(valX,2)        
        v = [valX(1,i) - valX(2,i); valY(1,i) - valY(2,i)];    % direction-vector
        w = [-v(2,1);v(1,1)] / (norm(v)*90);       % v*w = 0
        v = v / 4;                          
       
        % add or substract w and v to change the position of the line
        valX(1,i)           = valX(1,i) - w(1,1) - v(1,1);
        valX(2,i)           = valX(2,i) - w(1,1) + v(1,1);
    
        valY(1,i)           = valY(1,i) - w(2,1) - v(2,1);
        valY(2,i)           = valY(2,i) - w(2,1) + v(2,1);
    end
    
    % no color   
    valC = zeros(size(valX));
end

% function which returns the input parameters of the patch function to
% draw red bold lines for marked sides
function [valX, valY, valC] = getMarkXYC(c4n, n4sM)
    % coordinates for all nodes of marked sides
    X1 = c4n(n4sM(:,1), 1);
    X2 = c4n(n4sM(:,2), 1);
    Y1 = c4n(n4sM(:,1), 2);
    Y2 = c4n(n4sM(:,2), 2);   
    
    % combined coordinates in patch-style (would be lines up to here)
    X = [X1';X2'];
    Y = [Y1';Y2'];    
    valX = zeros(4,size(X,2));
    valY = zeros(4,size(Y,2));

    % makes boxes out of the edges so they can be seen better
    for i=1 : size(X,2)
        v = [X(1,i) - X(2,i); Y(1,i) - Y(2,i)];    % direction-vector
        w = [-v(2,1);v(1,1)] / (norm(v)*100);      % v*w = 0             
        
        % add or substract w to make it bold
        valX(1,i) = X(1,i) - w(1,1);        
        valX(2,i) = X(2,i) - w(1,1);
        valX(3,i) = X(2,i) + w(1,1);
        valX(4,i) = X(1,i) + w(1,1);
    
        valY(1,i) = Y(1,i) - w(2,1);
        valY(2,i) = Y(2,i) - w(2,1);
        valY(3,i) = Y(2,i) + w(2,1);
        valY(4,i) = Y(1,i) + w(2,1);
    end
    
    %color: red
    valC = zeros(1,size(valX,2),3);
    valC(1,:,1) = ones(1,size(valX,2),1);    
end

% function which returns the input parameters of the patch function to
% draw RGB-colored triangles (number of marked sides)
function [valX, valY, valC] = getColoredXYC(c4n, n4e, n4sR)    
    % coordinates for triangles 
    X1 = c4n(n4e(:,1), 1);
    X2 = c4n(n4e(:,2), 1);
    X3 = c4n(n4e(:,3), 1);
    Y1 = c4n(n4e(:,1), 2);
    Y2 = c4n(n4e(:,2), 2);
    Y3 = c4n(n4e(:,3), 2);    
    
    % combined coordinates in patch-style
    valX = [X1';X2';X3'];
    valY = [Y1';Y2';Y3'];
            
    % numbers of Nodes and Elements
    nrNodes = size(c4n,1);
    nrElems = size(n4e,1);
    valC = ones(1,nrElems,3);    %default color: white
    
    % compute newNodes4n to find new nodes faster as in refineRGB.m
    newNode4n = sparse(n4sR(:,1),n4sR(:,2),(1:size(n4sR,1))'+ nrNodes, nrNodes, nrNodes);
    newNode4n = newNode4n + newNode4n';
    
    for curElem = 1 : nrElems
        % the three nodes of the current element
        curNodes = n4e(curElem,:);
        curNewNodes = [newNode4n(curNodes(1),curNodes(2));
                       newNode4n(curNodes(2),curNodes(3));
                       newNode4n(curNodes(3),curNodes(1));
                      ];
        if nnz(curNewNodes) == 1     % green if one side is marked
            valC(1,curElem,:) = [0 1 0];        
        elseif nnz(curNewNodes) == 2 % blue if two sides are marked
            valC(1,curElem,:) = [0 0 1];   
        elseif nnz(curNewNodes) == 3 % red if all sides are marked
            valC(1,curElem,:) = [1 0 0];          
        end
    end
end

function waitForUser()
    input(sprintf('Enter to continue: '));
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
