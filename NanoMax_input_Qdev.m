% In order to open saved objects, the corresponding class definition has to
% be available

clear;
close all;

% add class definitions and functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('functions');
addpath('functions/openFunctions');
addpath('classDef');
addpath('gui');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAXWELL cluster / DESY
% input_param.pre_path          = fullfile('/data','netapp','dzhigd','Experiments','34_IDC_APS_STO_AFM_run3_2018','Dmitry_APS_34ID-C_20181106','data');
% input_param.save_folder       = fullfile('/data','netapp','dzhigd','Experiments','34_IDC_APS_STO_AFM_run3_2018','analysis','processed');
% input_param.white_field_path  = fullfile('/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/APS_34IDC_viewer','white_field_1.mat');
% input_param.dark_field_path   = fullfile('/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/APS_34IDC_viewer','dark_field_1.mat');
% input_param.mask_path         = '/data/netapp/dzhigd/Experiments/34_IDC_APS_STO_AFM_run3_2018/analysis/processed/Sample1_B_pristine_mask.mat';
% input_param.beamtime_id       = 'Dzhigaev1118';
% input_param.sample_name       = 'Sample1';

% work PC LUND / check the drive name! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All folder names start with capital letters
% File structure:
% 
% Data
% - DRIVE -> USER -> 'Data' -> beamtime_id -> [beamline dependent structure]
% 
% Analysis
% - DRIVE -> USER -> 'Projects' -> project_name -> 'Experiments' ->
% beamtime_id -> 'Analysis' -> 'Processed' -> location_id ['Onsite'/'Offsite'] 
% 
% Masks should be located in respective folders of the scan for single
% scans or in core analysis folder in case of multiple scans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select beamline to work-with %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 34IDC (APS)
% NanoMax (MAXIV)
% P10 (PETRAIII)
input_param.beamline_id       = 'NanoMax';

% Data Input Routing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_param.pre_path_data       = fullfile('E:','DDzhigaev','Data'); % General data folder
input_param.beamtime_id         = 'NanoMax_TFET_03_2019';
input_param.beamtime_prefix     = 'Dzhigaev1118';
input_param.sample_name         = 'Sample3';
input_param.scan_number         = 19;
input_param.detector_id         = 'merlin';
input_param.master_file_nanomax = 'C:\DDzhigaev\Projects\NanowireDevices\Experiments\NanoMax_TFET_03_2019\Analysis\Processed\Offsite\Sample_Qdev869_verticalMount\sample_Qdev896_horizontal.h5';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Misc Data Routing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_param.pre_path_analysis = fullfile('C:','DDzhigaev','Projects');
input_param.project_name      = 'NanowireDevices';
input_param.location_id       = 'Offsite'; % 'Onsite'
input_param.save_folder       = fullfile(input_param.pre_path_analysis,...
    input_param.project_name,'Experiments',input_param.beamtime_id,...
    'Analysis', 'Processed',input_param.location_id);
input_param.white_field_path  = fullfile(input_param.save_folder,'white_field_1.mat');
input_param.dark_field_path   = fullfile(input_param.save_folder,'dark_field_1.mat');
input_param.mask_path         = 'mask.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_param.exposure          = 1; % seconds <- has to come from SPEC like structure / fio-in case of P10
input_param.dead_time         = 0; % 2e-6 seconds
input_param.low_cutoff        = 0; % cutoff level of signal for general cleaning

mkdir(input_param.save_folder);

% - The object has to be created by opening method - %
scan = Scan(input_param);

% Reading
user_path = 'C:\DDzhigaev\Projects\NanowireDevices\Experiments\NanoMax_TFET_03_2019\Analysis\Processed\Offsite\Sample_Qdev869_verticalMount\';

scan.create_path(user_path);
% scan.read_tif;
scan.read_nanomax_merlin;

% scan.read_dark_field;
% scan.read_white_field;

% Correction bunch
% scan.correct_dark_field;
% scan.normalize_exposure; 
% scan.correct_white_field;
% scan.correct_low_cutoff;

% Mask stuff
% scan.create_mask;
% scan.save_mask;
scan.read_mask;
scan.correct_;

% Final hot/bad pixel corrections based on average
% data.correct_hot_pixel(100,30,1);
% data.correct_hot_pixel(91,106,1);
% data.correct_hot_pixel(123,58,1);
% data.correct_hot_pixel(1,18,1);

% Combining scans
% scan.combine;
% scan.crop;

scan.show_data_average('log'); % 'lin' - default, 'log'
scan.show_data_sum('log'); % 'lin' - default, 'log'

scan.save_data; % add user-name as an imput if you want a different name

disp('### Data import done! ###');

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
data.maximize;
data.minimize;
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

% Reconstruction -> plan OOP realisation
PHASOR_GUI;

% Visualisation/plotting
data.show_3d;
data.show_data_scroll('log'); % (<'lin','log'-default>, <max_value>)
data.show_data_single('log',1); % 'lin' - default, 'log'
data.show_data_average; % 'lin' - default, 'log'
data.show_data_sum; % 'lin' - default, 'log'
data.show_data_max;
data.show_data_min;
data.show_data_integral;
data.show_dark_field;
data.show_white_field;

% Output
data.save_gif;
data.save_mask;
data.save_data; % add user-name as an imput if you want a different name

close all;
disp('Test succesful!')
