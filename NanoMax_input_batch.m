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
% beamtime_id -> 'Analysis' -> 'Processed' -> location_id
% ['Onsite'/'Offsite'] -> sample_id
% 
% Masks should be located in respective folders of the scan for single
% scans or in core analysis folder in case of multiple scans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Batch initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_numbers = [102,103,106,114,115,117,118,119,121,123,128,130,134,135,136];

% Variable parameters go here
gonphis      = [16.2,16.18,16.12,15.98,15.96,15.92,15.9,15.88,15.84,15.8,15.7,15.66,15.58,15.56,15.54];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for scanN = 1:numel(scan_numbers)
    fprintf('Current scan processing: %d\n',scan_numbers(scanN));
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
    input_param.scan_number         = scan_numbers(scanN);
    input_param.detector_id         = 'merlin';
    input_param.master_file_nanomax = '/data/netapp/dzhigd/Experiments/NanoMax_TFET_03_2019/Analysis/Processed/Offsite/Sample_Qdev869_verticalMount/sample_Qdev896_horizontal.h5';
    input_param.save_diff           = 0;
    
    % input_param.fluo_path           = '/data/netapp/dzhigd/Experiments/NanoMax_TFET_03_2019/Analysis/Processed/Offsite/Sample_Qdev869_verticalMount/scan_0019_xspress3_0000.hdf5';
    
    % Crop settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_param.crop_flag         = 0;

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

    input_param.nanomax.gonphi          = gonphis(scanN); % [deg] - can be a range of angles
    input_param.nanomax.gontheta        = 0; % [deg] - can be a range of angles

    input_param.nanomax.radius          = 0.4; % [m]
    input_param.nanomax.delta           = 0.5092; % [deg] rotation of detector around horizontal axis of the sample
    input_param.nanomax.gamma           = 33; % [deg] rotation of detector around vertical axis of the sample
    input_param.nanomax.detector_pitch  = 55e-6; % [m]
    input_param.nanomax.wavelength      = 15000; % [eV]

    input_param.nanomax.direct_beam     = [125,275]; % [h,v]

    % Scan parameters
%     input_param.nanomax.step_h          = sind(input_param.nanomax.gonphi)*(63-45)*1e-6/85;
    input_param.nanomax.step_h          = (63-45)*1e-6/85;
    input_param.nanomax.step_v          = 5e-6/130;
    
    % Ring parameters
    input_param.nanomax.read_beam_current = 0;
    
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
        'Analysis', 'Processed',input_param.location_id,'Qdev');
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
    
    % Correction bunch


    % Mask stuff
    scan.read_mask;
    scan.correct_mask_nanomax;

    % Final hot/bad pixel corrections based on average


    % Combining scans

    scan.data = permute(scan.data,[2,1,3,4]);
    scan.flip('lr'); 
    scan.save_data; % add user-name as an imput if you want a different name
    
    disp('### Data import done! ###');
    clear scan;
end

%% Processing


%% Addons


%% Methods testing
% The full testing with all parameters was not performed yet
% 
% Reading
scan.read_tif; % can be used as a revert all 
scan.read_mask;
scan.read_dark_field;
scan.read_white_field;

% Processing
scan.average; % option is the final dimensionality of the output
scan.maximize;
scan.minimize;
scan.integrate;
scan.prepare_3d;
scan.create_mask; % options: '2D', '3D' - default
scan.combine; % works for 4D arrays
scan.flip('lr') % options: 'lr' - left-right, 'ud' - up-down
scan.bin2D(binning_size); % average binning_size number of pixels in the data

% Correction
scan.correct_dark_field;
scan.normalize_exposure; 
scan.correct_white_field;
scan.correct_low_cutoff;
scan.correct_mask;
scan.correct_mask_nanomax;
scan.correct_hot_pixel(100,30,0); % input: <[x,y] of the pixel, interpolate - 0/1>

% Reconstruction -> plan OOP realisation
PHASOR_GUI;

% Visualisation/plotting
scan.show_3d;
scan.show_data_scroll('log'); % (<'lin','log'-default>, <max_value>)
scan.show_data_single('log',1); % 'lin' - default, 'log'
scan.show_data_average; % 'lin' - default, 'log'
scan.show_data_sum; % 'lin' - default, 'log'
scan.show_data_max;
scan.show_data_min;
scan.show_data_integral;
scan.show_dark_field;
scan.show_white_field;

% Output
scan.save_gif;
scan.save_mask;
scan.save_data; % add user-name as an imput if you want a different name

close all;
disp('Test succesful!')
