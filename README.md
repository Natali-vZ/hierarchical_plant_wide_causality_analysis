# Hierarchical plant-wide causality analysis
Code for the plant-wide causality map in a  hierarchical approach for causality analysis, as described in:
(1) MEng thesis: van Zijl, N. (2020). Improving the interpretability of causlity maps for fault identification, Stellenobsch University.
https://scholar.sun.ac.za/bitstream/handle/10019.1/109269/vanzijl_improving_2020.pdf?sequence=1&isAllowed=y
(2) Journal Paper: van Zijl, N., Bradshaw, S. M., Auret, L., & Louw, T. M. (2021). A Hierarchical Approach to Improve the Interpretability 
of Causality Maps for Plant-Wide Fault Identification. Minerals, 11(8), 823. https://doi.org/10.3390/min11080823

For a summary, see the presentation, 'Final feedback_02-09-2020' in this repo.
This code follows the approach named 'PS-PC1' in the presentation, and
incorporates the following tools covered in the presentation:
(1) Incoporating process knowledge by validating data-based connections
with a connectivity matrix.
(2) Incorporating process knowledge by constraining potential root causes.
(3) Tools for interpretation: Display node rankings

This code makes use of publically available toolboxes and repositories:
(1) The Multivariate Granger Causality (MVGC) Toolbox: https://uk.mathworks.com/matlabcentral/fileexchange/78727-the-multivariate-granger-causality-mvgc-toolbox
(2) IFAC-LoopRank: https://github.com/ProcessMonitoringStellenboschUniversity/IFAC-LoopRank

The main script is "MAIN_Hierarchical", and "Sim_PlantModel_repIDs", "PLANTWIDE_ALL_DATA_SmallerFreqOsc_SmallerExogDist3", and 
"PLANTWIDE_Connectivity_matrix_FINAL" contain the simulated case study data and process knowledge.
