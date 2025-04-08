% Copyright (C) 2020 N van Zijl

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

% Improving the interpretability of causality maps for fault identification
% MEng thesis, Stellenbosch University

% N van Zijl
% Department of Process Engineering
% University of Stellenbosch

%%
% *Note 1: Connectivity matrix should consist of equipment (units & pipelines),
% controllers, sensors & valves.

% *Note 2: Set up for a causality matrix consisting of 1 PC from each section (where each section has an equal number of variables). 
% Adapting to allow for groups(sections) to have different nummber of
% variables

%% Note from code for CL causality matrix: 
% IDs for causality
% matrix will only be controller names, which have to correspond to
% controller names in connectivity matrix. Part of code the code that will
% have to be adapted for normal causality matrix, is marked with '***' in a
% comment above and below it.

% *NOTE: The way this fcn is now, it will only work for i PC from each
% section. Because numSections is determined as a dimension of
% causalM_original. If you want to update it to allow for more PCs per
% section, then you'd have to specify the number of sections earlier on and
% feed that as an input to this fcn.

function causalM_refined = validate_propPath(causalM_original, connectM, all_IDs_causal, IDs_connect)
    [numGroups, ~] = size(causalM_original);
    numVars = NaN(1, numGroups); % Row vector that contains the no. of vars in each group
    for i = 1:numGroups
        [~, numVars(i)] = size(all_IDs_causal{1, i}); % Okay for now, because there are an equal number of variables in each section
    end
    causalM_refined = zeros(numGroups);
    % Step through causal matrix to check each connection:
    for i_sourceSection = 1:numGroups
        for i_sinkSection = 1:numGroups
            if i_sourceSection == i_sinkSection % exclude diagnoal to omit self loops (self loops are terrible for interpretability)
              causalM_refined(i_sourceSection, i_sinkSection) = 0;
            else 
                % ***
                % Find indices of entries where there are causal connections:
                if causalM_original(i_sourceSection, i_sinkSection) ~= 0
                    
                    % ***
                    % Instead of stepping through all the variables in the
                    % source section, just include all those variables in
                    % the list of start nodes. sn has to be a scalar, so
                    % will need to step throuhg all the variables after
                    % all:
                    for i_sourceSectionVars = 1:numVars(i_sourceSection)
                    ID_sn = all_IDs_causal{1, i_sourceSection}{i_sourceSectionVars}; % numVars must equal the number of variables in the current source section; generic for now, because all the sections have an equal number of variables
                    idx = find(strcmp(IDs_connect, ID_sn));
                    % ***

                    % What is the instrument connected to? (Check row & col
                    % corresponding to this instrument for 1's):
                    idx_row = find(connectM(idx, :)==1); % Row containing indices where there is a 1
                    idx_col = find(connectM(:, idx)==1); % Col containing indicies where there is a 1

                    % Save indices as start nodes (doesn't matter if they were 
                    % row or col indices, because variables are in the same order 
                    % in the rows & cols, so can use the obtained idices directly 
                    % as row indices for the start nodes):
                    snVect = zeros(size(idx_row, 2)+size(idx_col, 1), 1); %snVect is a col vector
                    snVect(1:size(idx_row, 2)) = idx_row;
                    snVect(size(idx_row, 2)+1:size(idx_row, 2)+size(idx_col, 1)) = idx_col;

                    % DFS:
                    % Create digraph out of connectivity matrix
                    connectM_dg = digraph(connectM, IDs_connect); % dg = digraph

                    % Pre-allocate space to store the paths
                    connectM_paths = cell(1, size(snVect, 1));

                    % DFS to get all nodes that can be reached from the specified
                    % start nodes
                    for i_DFS = 1:size(snVect, 1)
                       connectM_paths{1, i_DFS} = dfsearch(connectM_dg, IDs_connect{1, snVect(i_DFS)}); % 1st entry contains nodes reachable from 1st specified start node, and so on..
                    end

                    % Is the target node in the list of nodes that can be
                    % reached from the start nodes? If yes, causal connection
                    % is validated. If no, causal connection is deemed spurious
                    % and removed.
                    
                    % IDs_causal(i_causal_col) is the target node -
                    % Make it all nodes in the sink section
                        
                    %The for loop for the variables in the sink section
                    %must be here:
                    b = 0;
                    for i_sinkSectionVars = 1:numVars(i_sinkSection)
                        for i_targetReachable = 1:size(connectM_paths, 2)
                            a = strcmp(all_IDs_causal{1, i_sinkSection}(i_sinkSectionVars), connectM_paths{1, i_targetReachable});
                            b = b+sum(a);
                            if b > 0
                              causalM_refined(i_sourceSection, i_sinkSection) = causalM_original(i_sourceSection, i_sinkSection);
                              break;
                            end  
                        end 
                    end 
                    end 
                end
            end
        end 
    end
end