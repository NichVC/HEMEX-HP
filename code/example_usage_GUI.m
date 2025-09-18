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
load ../data/Pig_kidney_demo.mat

%% Call the HEMEX GUI function and proceed to following the guide in "docs"
% Standard model parameter bounds are defined in the function, and can thus
% be omitted. However, it can also be given as an additional parameter.
HEMEX_GUI_main(data, flip_P, flip_L, TR);

% The bounds can be given in this form:
% lower bound | start guess | upper bound
bounds = [0 0.01 1;... % kpl (pyrurvate-to-lactate conversion rate)
          0 0 0;... % klp (lactate-to-pyruvate conversion rate)
          0.01 1/30 0.05;... % rp (1/T1 pyruvate relaxation)
          0.8 1 1.2;... % rl (1/T1 lactate relaxation) [scaled from rp]
          0.05 1 20;... % k (pyruvate permability) [scaled from kpl]
          0 0 10;... % t0_delay (time-delay between AIF and voxel data)
          0.5 5 30;... % mu (input parameter to the hemodynamic residue function)
          0.5 3 30]; % sigma (input parameter to the hemodynamic residue function)
% The function is then called by:
% HEMEX_GUI_main(data, flip_P, flip_L, TR, bounds);