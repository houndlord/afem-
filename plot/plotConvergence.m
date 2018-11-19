function plotConvergence(nrDoF4lvl, error4lvl, OPTname)
%% Plot error for different levels.
%   plotConvergence(nrDoF4lvl, error4lvl) plots the error given by
%           error4lvl over the degrees of freedom given by nrDoF4lvl, both
%           using logarithmic scale. The optional String argument together
%           with the convergence rate is added to the legend.
    if nargin > 2
        name = OPTname;
    else
        name = '';
    end

    cOrder = [0 0 1; 1 0 0; 0 .7 0; .7 0 .7; .7 .7 0; 0 0 0];

    nlines = get(gcf, 'UserData');

    if ~isempty(nlines)
      nlines = nlines + 1;
    else
      nlines = 1;
    end
    set(gcf, 'UserData', nlines);

    %% Plot.
    plot = loglog(nrDoF4lvl, error4lvl, '-s', 'Color', cOrder(nlines, :));
    title('Convergence history plot');

    %% Set name and title of the figure.
    set(plot,'DisplayName',name);

    legend('off');
    legend('show');
    handle = findobj(gcf,'type','axes','Tag','legend');
    set(handle,'Location','best');

    drawnow;
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
