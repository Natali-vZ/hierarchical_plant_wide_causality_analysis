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
% MEng thesis, Stellenbosch University &
% van Zijl, N., Bradshaw, S. M., Auret, L., & Louw, T. M. (2021). 
% A Hierarchical Approach to Improve the Interpretability of Causality Maps 
% for Plant-Wide Fault Identification. Minerals, 11(8), 823. 
% https://doi.org/10.3390/min11080823

% N van Zijl
% Department of Process Engineering
% University of Stellenbosch

%% Description:
% Code for the plant-wide causality map in a  hierarchical approach for
% causality analysis.
% For more info, see the presentation, 'Final feedback_02-09-2020':
% This code follows the approach named 'PS-PC1' in the presentation, and
% incorporates the following tools covered in the presentation:
% (1) Incoporating process knowledge by validating data-based connections
% with a connectivity matrix.
% (2) Incorporating process knowledge by constraining potential root
% causes.
% (3) Tools for interpretation: Display node rankings

%% General notes for using this code:
% (1) Sections of code that you need to adapt for your application say 'EDIT'
% in their headings.
% (2) Sections of code that will not need to be adapted for every
% application, but that you may want to play around with, say 'POSSIBLY
% EDIT' in their headings.
% (3) Each section of code that needs/may need adaptation, contains notes
% right at the beginning of the section, explaining how to edit the code.
% (4) When entire sections of code need to be replaced with your own code
% to adapt to your application, I marked my code with '% ~' at the 
% beginning and end of the part that needs to be replaced.


%% Initialise
clc; clear; close all; 

%% EDIT: Load data
% Notes: 
% (1) Load the data in a cell array named 'ALL_DATA'.
% (2) In the matrices inside the cell array, rows are observations & columns are variables.
% (3) Make sure to sub-sample if necessary.

% ~ 
load 'PLANTWIDE_ALL_DATA_SmallerFreqOsc_SmallerExogDist3.mat' 
data_1 = [data_1(:, 1:4) data_1(:, 7)]; 
data_2 = [data_2(:, 1:4) data_2(:, 7)];
data_3 = [data_3(:, 1:4) data_3(:, 7)];
% Cut data to remove "start-up" period in my sim:
cut_data = 1800;
% Sample from data
ST = 5;
TW = 5*602;
Data_1 = data_1(cut_data:ST:TW, :); 
Data_2 = data_2(cut_data:ST:TW, :);
Data_3 = data_3(cut_data:ST:TW, :);
ALL_DATA = {Data_1, Data_2, Data_3};
% ~

%% EDIT: Load IDs 
% Notes:
% (1) The IDs are loaded as row string vectors (each representing a plant section)
% in a cell array named 'all_IDs_forVars'.
% (2) Make sure that the order of these IDs corresponds to the order of the
% variables in the data loaded above.
% (3) Make sure that these IDs correspond to the tag names present in the connectivity matrix.

% ~
IDs_1 = {'LI-001', 'LI-002', 'FI-001', 'FI-002', 'FI-003'};
IDs_2 = {'LI-003', 'LI-004', 'FI-004', 'FI-005', 'FI-006'};
IDs_3 = {'LI-005', 'LI-006', 'FI-007', 'FI-008', 'FI-009'};
all_IDs_forVars = {IDs_1, IDs_2, IDs_3}; % IDs grouped according to plant sections
% ~

%% EDIT: Load connectivity matrix & its IDs:
% Notes:
% (1) The connectivity matrix is loaded WITHOUT lables as a square matrix named 'connectM'.
% (2) The connectivity matrix IDs are loaded as a row string vector named
% 'IDs_connect'.

% ~
connectM = xlsread('PLANTWIDE_Connectivity_matrix_FINAL.xlsx','Sheet1', 'B2:BK63');
[~, IDs_connect, ~]= xlsread('PLANTWIDE_Connectivity_matrix_FINAL.xlsx','Sheet1','B1:BK1'); % IDs for everything in connectM
% ~

%% POSSIBLY EDIT: Parameters for GCCA mode in MVGC toolbox:
% Notes: 
% (1) This can run as it is now, but you can play with the parameters if
% you want to.
% (2) You will likely only edit 'momax', which is the maximum allowable
% model order for the AR models, and 's.alpha', which is the significance
% level for the statistical test.
% (3) Typical values: 
% momax: I honestly think 60 is a good max for plant-wide causality
% analysis.
% s.alpha: 0.01 or 0.05 are typical values. You can also set a very strict
% threshold by using 0 (which equates to 1e-16 in MATLAB).

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  
momax     = 60;
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
tstat     = 'F'; 

% Initial (default) values for sliders:
s.alpha     = 0.01; % or 0.05 or 0 (which equates to 1e-16)

%% POSSIBLY EDIT: To perform PCA:
% Notes: 
% (1) You might have to edit 'allIDs_forReps'. It is a matrix/vector that
% provides the IDs (according to groups/plant sections) on the plant-wide
% causality map. As it is, you can have up to 8 groups/plant sections.

numGroups = size(ALL_DATA, 2);
numPCs = 1; % At the moment this does nothing, because the 'pca_on_each_section' fcn has not been updated so that you can spcify the number of PCs to include; it includes 1 PC (unless you manually alter the fcn) 
level1_data = pca_on_each_section(ALL_DATA, numPCs);
X = level1_data';
x = X;
[nvars, nobs, ntrials] = size(x);

% IDs for PCs:
load 'Sim_PlantModel_repIDs.mat' % Can use this for up to 8 plant sections.
allIDs_forReps = allIDs_forPCs;

%% Model order & GC (for default settings):
p = estimateModelOrder(x, momax, icregmode, morder);
[F, sig, pval] = gc(x, s.alpha, p, regmode, mhtc);

%% Initial causality matrix with default alpha (0.05) & determine reachbility matrix:
causalM_sig = causalM_significant(F, sig);

% Reachability matrix:
R = reachabilityM(sig'); % Feeding transpose of sig, because sig is an output from the GCCA toolbox, so rows are variables & cols observations, but all my code is the other way around

%% Hybrid approach: Check CA results with connectivity matrix - Show connections not corroborated by connectivity matrix in grey (faded)
causalM_refined = validate_propPath(causalM_sig, connectM, all_IDs_forVars, IDs_connect); % This is where I'll have to check that the correct ID names are being used!

%% Plot initial causality map:
g = figure;
g = plot_fig(causalM_sig, allIDs_forReps(1:numGroups, 1), g, causalM_refined, R); % Here I specify that only the IDs for the 1st PCs/KPIs must be used/shown on the causlaity map

%% Causality matrix showing only statistically significant connections:
function [causalM_sig] = causalM_significant(F, sig)
causalM_sig = zeros(size(F));
[rows, cols] = size(causalM_sig);
for source_i = 1:rows
    for sink_i = 1: cols
        if sig(sink_i, source_i) == 1
            causalM_sig(source_i, sink_i) = F(sink_i, source_i);
        end 
    end 
end 
end

%% Determine reachability matrix:
function R = reachabilityM(sig)
[nvars, ~] = size(sig);
adjM = sig; % Adjacency matrix
adjM(isnan(adjM))=0;

R = zeros(nvars); % Preallocate space for reachability matrix
adjM_power = cell(nvars, 1); % Preallcoate space for powers of adjacency matrix
adjM_power{1} = adjM; % Set first element of cell containing power of adjacency matrix equal to adjacency matrix
for i = 2:nvars
    adjM_power{i} = adjM^i + adjM_power{i-1};
end
for i_R_rows = 1:nvars
    for i_R_cols = 1:nvars
        if adjM_power{nvars, 1}(i_R_rows, i_R_cols) ~= 0
            R(i_R_rows, i_R_cols) = 1;
        else
            R(i_R_rows, i_R_cols) = 0;
        end
    end
end 
end

%% Model order estimation - 'MVGC toolbox' (same for GCCA & recommended pathway)
function bmo_AIC = estimateModelOrder(X, momax, icregmode, morder)

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

[~,bmo_AIC] = min(AIC);
[~,bmo_BIC] = min(BIC);

% Plot information criteria.
iac_plot = figure;
plot((1:momax)',[AIC BIC]);
legend('AIC','BIC');

fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
fprintf('best model order (BIC) = %d\n',bmo_BIC);

% Select model order
if strcmpi(morder,'AIC')
    morder = bmo_AIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
else
    morder = bmo_BIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
end
end

%% Granger causality estimation - GCCA mode in 'MVGC toolbox'
function [F, sig, pval] = gc(X, alpha, morder, regmode, mhtc) 
tstat = 'F';

[nvars, nobs, ntrials] = size(X);
ptic('\n*** GCCA_tsdata_to_pwcgc... ');
[F,A,SIG] = GCCA_tsdata_to_pwcgc(X,morder,regmode); % use same model order for reduced as for full regressions
ptoc;

% Check for failed (full) regression
assert(~isbad(A),'VAR estimation failed');

% Check for failed GC calculation
assert(~isbad(F,false),'GC calculation failed');

% Check VAR parameters (but don't bail out on error - GCCA mode is quite forgiving!)
rho = var_specrad(A);
fprintf('\nspectral radius = %f\n',rho);
if rho >= 1,       fprintf(2,'WARNING: unstable VAR (unit root)\n'); end
if ~isposdef(SIG), fprintf(2,'WARNING: residuals covariance matrix not positive-definite\n'); end

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.
pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat);
sig  = significance(pval,alpha,mhtc);
end

%% Plot causality map
function [g] = plot_fig(causalM_sig, IDs, g, causalM_refined, R)
G = digraph(causalM_sig, IDs, 'omitselfloops');
figure(g);
GPlot = plot(G, 'Layout', 'layered', 'MarkerSize', 8);

% Set edge colour to black:
GPlot.EdgeColor = [0 0 0];
% Make arrowheads bigger & place them closer to the sink nodes:
GPlot.ArrowSize = 15;
GPlot.ArrowPosition = 0.9; % [0 1], where 0 is right next to source, and 1 is right next to sink

% Colour nodes according to NodeRank:
ranking_John_G = LoopRank(causalM_refined);
G.Nodes.NodeColors = ranking_John_G;
GPlot.NodeCData = G.Nodes.NodeColors;
colormap jet
colorbar

if nargin > 4    
    % Fade nodes that cannot reach all variables in specified group
    vars_toReach = IDs; % For now, the group of variables that a node needs to reach to not be faded, includes all variables.
    % If there are no connections to a node, it should be excluded from
    % this group that must be reachable - surely?
    adj = adjacency(G);
    for i_NoEdges = 1:size(vars_toReach, 1)
        if sum(adj(i_NoEdges, :))==0 && sum(adj(:, i_NoEdges))==0
            vars_toReach{i_NoEdges} = 'delete';
        end 
    end 
    vars_toReach = vars_toReach(~strcmp(vars_toReach, 'delete'));
    
    nvars_toReach = size(vars_toReach, 1)-1; % Minus one, because the variable does not need to reach itself.1
    [nvars, ~] = size(causalM_sig);
    R_sum = zeros(nvars, 1); % Pre-allocate space; difficult to explain what this is for - just check the code :)
    for i_source = 1: nvars
        for i_sink = 1: nvars
            if i_source ~= i_sink % Exclude diagonal
                R_sum(i_source) = R_sum(i_source) + R(i_source, i_sink);
            end
        end
    end
    idx_to_fade = find(R_sum<nvars_toReach);
    nodes_to_fade = IDs(idx_to_fade);
    
    % Plot a subgraph (containing the faded nodes) over the graph:
    fadedG = subgraph(G, nodes_to_fade);
    idx_edgesRemove = 1:numedges(fadedG); % Remove all edges from subgraph, so that when it gets plotted over the original graph, edges don;t show up twice.
    fadedG_NodesOnly = rmedge(fadedG, idx_edgesRemove);
    hold on
    fadedGPlot = plot(fadedG_NodesOnly, 'XData', GPlot.XData(idx_to_fade), 'YData', GPlot.YData(idx_to_fade), 'MarkerSize', 8, 'NodeColor', [227/255 225/255 226/255]);
    % Remove node labels (so the original ones will show):
    fadedGPlot.NodeLabel = {};
end     

% Make node labels (of original graph) bigger:
labelSizes = 20*ones(1, length(IDs));
GPlot.NodeLabel = {};
hold on
for i = 1:length(IDs)
    text(GPlot.XData(i), GPlot.YData(i), IDs(i), 'FontSize', labelSizes(i));
end
hold off

if nargin > 3 % If causalM_refined is also fed as input (i.e. if process knowledge is to be included in the causality map)
    % Obtain set of edges that are not validated by connectivity matrix:
    causalM_notVal = causalM_sig - causalM_refined;
    causalM_notVal_dg = digraph(causalM_notVal, IDs);
    edgesNotVal = causalM_notVal_dg.Edges.EndNodes;
    s = edgesNotVal(:, 1);
    t = edgesNotVal(:, 2);
    s = s';
    t = t';
    % Display the edges that are not validated by the connectivity matrix with dotted lines:
     highlight(GPlot, s, t, 'LineStyle', '--');
end

% Set edge thickness according to edge weights
if numedges(G) ~= 0
    G.Edges.LWidths = 7*G.Edges.Weight/max(G.Edges.Weight); % Scale line widths s.t. edge with greatest weight has line width of 7 (because edge line width greater than 7 starts becoming cumbersome).
    GPlot.LineWidth = G.Edges.LWidths;    
end
set(gcf, 'color', 'w');
end

%% PCA
function level1_data = pca_on_each_section(ALL_DATA, numPCs) % Can add functionality to spefify the number of PCs to include later.
% Only keep non-empty groups:
ALL_DATA = ALL_DATA(~cellfun('isempty',ALL_DATA));

numGroups = size(ALL_DATA, 2);

[numObs, numVars] = size(ALL_DATA{1, 1}); % When the sections contain different amount of variables, I'll have to make a for loop & cell here as well, to get the number of variables in each section.

ALL_DATA_z = cell(size(ALL_DATA)); % Cell to store all standardised datasets in
coeff_score = cell(2, size(ALL_DATA, 2)); % Cell to store all coefficients & scores from PCA for each dataset in
for i = 1:numGroups
    % Standardise data:
    ALL_DATA_z{1, i} = zscore(ALL_DATA{1, i});
    % Perform PCA on each module/plant section matrix:
    [coeff, score, ~] = pca(ALL_DATA_z{1, i});
    coeff_score{1, i} = coeff;
    coeff_score{2, i} = score;
end
% Set up matrix with first PC from each module/plant section:
level1_data = NaN(numObs, numGroups);
for i = 1:numGroups
    level1_data(:, i) = coeff_score{2, i}(:, 1);
end
end
