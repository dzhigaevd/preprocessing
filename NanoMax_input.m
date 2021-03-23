% In order to open saved objects, the corresponding class definition has to
% be available

clear;
close all;

addpath('openFunctions');
addpath('functions');
addpath('classDef');
addpath('gui');

input_param.pre_path          = fullfile('/data','netapp','dzhigd','Experiments','34_IDC_APS_STO_AFM_run3_2018','Dmitry_APS_34ID-C_20181106','data');
input_param.save_folder       = fullfile('/data','netapp','dzhigd','Experiments','34_IDC_APS_STO_AFM_run3_2018','analysis','processed');
input_param.white_field_path  = fullfile('/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/APS_34IDC_viewer','white_field_1.mat');
input_param.dark_field_path   = fullfile('/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/APS_34IDC_viewer','dark_field_1.mat');
input_param.mask_path         = '/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/analysis/processed/Sample1_B_pristine_mask.mat';
input_param.beamtime_id       = 'Dzhigaev1118';
input_param.sample_name       = 'Sample1';
input_param.scan_number       = [151,155,158,161,164,167];
input_param.maindirswitch     = 4;
input_param.exposure          = 1; % seconds <- has to come from SPEC like structure / fio-in case of P10
input_param.dead_time         = 0; % 2e-6 seconds
input_param.low_cutoff        = 1; % cutoff level of signal for general cleaning

% input_param.file_names        = createPath_aps(input_param);

mkdir(input_param.save_folder);

% - The object has to be created by opening method - %
data = Scan(input_param);

% Reading
data.create_path;
data.read_tif;
data.read_mask;
data.read_dark_field;
data.read_white_field;

% Correction bunch
data.correct_dark_field;
data.normalize_exposure; 
data.correct_white_field;
data.correct_low_cutoff;
% data.create_mask;
data.correct_mask;

% Final hot/bad pixel corrections based on average
data.correct_hot_pixel(100,30,1);
data.correct_hot_pixel(91,106,1);
data.correct_hot_pixel(123,58,1);
data.correct_hot_pixel(1,18,1);

% Combining scans


data.show_data_average('log'); % 'lin' - default, 'log'

%% Methods testing
% The full testing with all parameters was not performed yet

% Reading
data.read_tif; % can be used as a revert all 
data.read_mask;
data.read_dark_field;
data.read_white_field;

% Processing
data.average; % option is the final dimensionality of the output
data.sum;
data.integrate;
data.prepare_3d;
data.create_mask; % options: '2D', '3D' - default
data.combine; % works for 4D arrays
data.flip('lr') % options: 'lr' - left-right, 'ud' - up-down

% Correction
data.correct_dark_field;
data.normalize_exposure; 
data.correct_white_field;
data.correct_low_cutoff;
data.correct_mask;
data.correct_hot_pixel(100,30,0); % input: <[x,y] of the pixel, interpolate - 0/1>

% Reconstruction

% Visualisation/plotting
data.show_3d;
data.show_data_scroll('lin'); % (<'lin','log'-default>, <max_value>)
data.show_data_single('lin',40); % 'lin' - default, 'log'
data.show_data_average; % 'lin' - default, 'log'
data.show_data_sum; % 'lin' - default, 'log'
data.show_data_integral;
data.show_dark_field;
data.show_white_field;

% Output
data.save_gif;
data.save_mask;
data.save_data; % add user-name as an imput if you want a different name

close all;
disp('Test succesful!')
