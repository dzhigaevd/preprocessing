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
input_param.beamline_id       = 'nanomax';

% Data Input Routing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_param.pre_path_data       = fullfile('/data','netapp','dzhigd'); % General data folder
input_param.beamtime_id         = 'NanoMax_TFET_03_2019';
input_param.beamtime_prefix     = 'Dzhigaev1118';
input_param.sample_name         = 'Sample3';
input_param.scan_number         = 19;
input_param.detector_id         = 'merlin';
input_param.master_file_nanomax = 'C:\DDzhigaev\Projects\NanowireDevices\Experiments\NanoMax_TFET_03_2019\Analysis\Processed\Offsite\Sample_Qdev869_verticalMount\sample_Qdev896_horizontal.h5';

input_param.fluo_path           = '/data/netapp/dzhigd/Experiments/NanoMax_TFET_03_2019/Analysis/Processed/Offsite/Sample_Qdev869_verticalMount/scan_0019_xspress3_0000.hdf5';

% Crop settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_param.crop_flag         = 1;

input_param.start_row         = 39;
input_param.end_row           = 83;
input_param.start_column      = 21;
input_param.end_column        = 194;

input_param.roi               = [1,inf,1,inf]; % [roiYstart, roiYend, roiXstart, roiXend]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NanoMax convention:
% gamma - horizontal detector
% delta - vertical detector
% gonphi - rotation about vertical axis
% gontheta - rotation about horizontal axis

input_param.nanomax.gonphi          = 15.67; % [deg] - can be a range of angles
input_param.nanomax.gontheta        = 0; % [deg] - can be a range of angles

input_param.nanomax.radius          = 0.4; % [m]
input_param.nanomax.delta           = 0.5092; % [deg] rotation of detector around horizontal axis of the sample
input_param.nanomax.gamma           = 33; % [deg] rotation of detector around vertical axis of the sample
input_param.nanomax.detector_pitch  = 55e-6; % [m]
input_param.nanomax.wavelength      = 15000; % [eV]

input_param.nanomax.direct_beam     = [125,275]; % [h,v]

% Scan parameters
input_param.nanomax.step_h          = 83e-9;
input_param.nanomax.step_v          = sind(input_param.nanomax.gonphi)*0.2e-6;

% 34IDC convention:
% gamma - horizontal detector
% delta - vertical detector
% gonphi - rotation about vertical axis
% gontheta - rotation about horizontal axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Misc Data Routing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_param.pre_path_analysis = fullfile('/data','netapp','dzhigd');
input_param.project_name      = '';
input_param.location_id       = 'Offsite'; % 'Onsite'
input_param.save_folder       = fullfile(input_param.pre_path_analysis,...
    input_param.project_name,'Experiments',input_param.beamtime_id,...
    'Analysis', 'Processed',input_param.location_id);
input_param.white_field_path  = fullfile(input_param.save_folder,'white_field_1.mat');
input_param.dark_field_path   = fullfile(input_param.save_folder,'dark_field_1.mat');
input_param.mask_path         = '/data/netapp/dzhigd/Experiments/NanoMax_TFET_03_2019/Analysis/Processed/Offsite/mask.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_param.save_formats      = 'mat';

input_param.exposure          = 1; % seconds <- has to come from SPEC like structure / fio-in case of P10
input_param.dead_time         = 0; % 2e-6 seconds
input_param.low_cutoff        = 0; % cutoff level of signal for general cleaning

mkdir(input_param.save_folder);

% - The object has to be created by opening method - %
scan = Scan(input_param);

% Reading
user_path = '/data/netapp/dzhigd/Experiments/NanoMax_TFET_03_2019/Analysis/Processed/Offsite/Sample_Qdev869_verticalMount';

scan.create_path(user_path);
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
scan.correct_mask_nanomax;

% Final hot/bad pixel corrections based on average
% data.correct_hot_pixel(100,30,1);
% data.correct_hot_pixel(91,106,1);
% data.correct_hot_pixel(123,58,1);
% data.correct_hot_pixel(1,18,1);

% Combining scans
% scan.combine;
% scan.crop;

% scan.show_data_average('log'); % 'lin' - default, 'log'
% scan.show_data_sum('log'); % 'lin' - default, 'log'

scan.data = permute(scan.data,[2,1,3,4]);
scan.flip('lr'); 
scan.save_data; % add user-name as an imput if you want a different name

disp('### Data import done! ###');

%% Processing
load('/data/netapp/dzhigd/Experiments/NanoMax_TFET_03_2019/Analysis/Processed/Offsite/Sample3_19/Sample3_19.mat')
% scan.data = permute(scan.data,[2,1,3,4]);
% scan.flip('lr'); % options: 'lr' - left-right, 'ud' - up-down
scan.average([3,4]);
scan.integrate([1,2]);
scan.show_data_average;
scan.show_data_integral;

alpha_map = scan.data_integral'./max(scan.data_integral(:));
figure; aa = imagesc(scan.data_meta.nanomax.,vVector,alpha_map);axis xy
axis image; colormap bone;
xlabel('Scan position, [um]');ylabel('Scan position, [um]');


scan.showScanCOM();
scan.showScanSTXMLive;
scan.save_data; % add user-name as an imput if you want a different name

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
scan.correct_mask_nanomax;
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
