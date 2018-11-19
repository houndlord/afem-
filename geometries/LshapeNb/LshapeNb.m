function LshapeNb()
% LSHAPENB The L-Shape geometry with partial Neumann boundary. In the 
% following scheme 'D' depicts Dirichlet boundary whereas 'N' depicts 
% Neumann boundary.
%
%     NNNNNNNNNNNNN
%     N           N
%     N           N
%     N     DDDDDDD
%     N     D
%     N     D       
%     NNNNNND
%
% Example: 
% [c4n n4e n4sDb n4sNb] = loadGeometry('LshapeNb',1);
%
% See also LSHAPE, LSHAPEROT
