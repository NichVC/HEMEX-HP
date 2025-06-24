% @author: Nichlas Vous Christensen
% @email: nvc@clin.au.dk
% @phone: +45 23464522
% @organization: Aarhus University, The MR Research Centre
% June 2025

%% First add the entire repository to path
addpath(genpath('../'));

%% Load the supplied data example of dynamic 13C imaging of pig kidneys
% Includes:
% data: 5D matrix with dimensions (40 x 40 x 2 x 20 x 2) [x, y, z, time, metabolite] with metabolite = 1 being pyruvate
% flip_P = 8 degrees
% flip_L = 70 degrees
% TR: repetition time of 3 seconds
load ../data/Pig_kidney_example.mat

%% Call the HEMEX GUI function and proceed to following the guide in "docs"
% Standard model parameter bounds are defined in the function, and can thus
% be omitted. However, it can also be given as an additional parameter.
HEMEX_GUI_main(data,flip_P,flip_L,TR);

