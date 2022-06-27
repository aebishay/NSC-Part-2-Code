global n_seeds Finaltimestep Doses chemo_sensitivity CONVERT2MICRON coh_limit chmtx_limit d_w cancer_size  cancer_center ...
    has_cancer modelNum modelType alpha4chmtx beta4dist d_g FolderName1 FolderName2 

%%% Standard values
n_seeds =           [300 600 900 1200 1500 1800];               %Total number of paths generated
%Finaltimestep = [14999 24999 34999 49999 59999 85000]; %Num of steps for each cell in doses
Finaltimestep = [2999 4999 6999 9999 11999 17000];
Doses = 6;
% TIMESTEP =          ???;                %Experimentally calculated tick per real time
CONVERT2MICRON =    54;                   %Avg μm per pixel in our figure
coh_limit =         0.4;                  %coherency threshold
chmtx_limit =       0.1;                  %chemotaxis threshold
d_w =               0.2;                %Reference step size on white matter
cancer_size =       [20, 20, 20];         %radius [x y z]
chemo_sensitivity = 0;

%%% Injection Center 
modelType =         'Intranasal'; 
% modelType =         'Intracerebral'; 


%%% Cancer Location
% cancer_center =     [230, 300, 175];    % right, front putamen 
% cancer_center =     [120, 400, 150];    % right, back putamen 
cancer_center = [79 84 78];      % right, front putamen 


%%% Conditional
has_cancer =        1;
modelNum =          2;


%%% Variables
alpha4chmtx =       1;                  %alpha(beta=1)  parameter for WM vs chemo
beta4dist =         1;                  %beta(alpha=1)  parameter for distance
d_g =              0.2;                  %reference step size on grey matter


FolderName1 = '/Figures/CancerNA/';     %Save plots to this folder when there is no cancer (Must pre-exist)
FolderName2 = '/Figures/Cancer/';       %Save plots to this folder when there is cancer (Must pre-exist)