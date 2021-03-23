%% Importing scans

clear;
close all;
clc;

user_path = '/home/jlastow/Documents/';
data_path = '/lunarc/nobackup/users/dzhigd/Data/MAXIV/NanoMax/2019032008/raw/sample_GaSb_W';
nanomax_preprocessing_path = '/home/jlastow/Documents/nanomax_preprocessing';

% add class definitions and functions 
addpath(fullfile(nanomax_preprocessing_path,'functions'));
addpath(fullfile(nanomax_preprocessing_path,'functions/openFunctions'));
addpath(fullfile(nanomax_preprocessing_path,'classDef'));
addpath(fullfile(nanomax_preprocessing_path,'gui'));

% Batch initialization 
scan_numbers = 185:215;

% Variable parameters go here
gonphis      = [5 5.02 5.04 5.06 5.08 5.10 5.12 5.14 5.16 5.18 5.20 5.22 5.24 5.26 5.28...
                5.30 5.32 5.34 5.36 5.38 5.40 5.42 5.44 5.46 5.48 5.50 5.52 5.54 5.56...
                5.58 5.60]; 
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for scanN = 1:numel(scan_numbers)
    fprintf('Current scan processing: %d\n',scan_numbers(scanN));

    input_param.beamline_id         = 'nanomax';

    % Data Input Routing 
    input_param.pre_path_data       = data_path; % General data folder
    input_param.sample_name         = 'Sample1';
    input_param.scan_number         = scan_numbers(scanN);
    input_param.master_file_nanomax = fullfile(data_path, 'sample_GaSb_W.h5');
    input_param.save_diff           = 0;
    
    input_param.mask_path           = fullfile(user_path, 'data/mask.mat');
    

    % Experiment parameters 
    input_param.nanomax.gonphi          = gonphis(scanN);
    input_param.scan.type               = 'flyscan' %flyscan or dmesh
    input_param.crop_flag               = 0;
   
    % Misc Data Routing 
    input_param.save_folder       = fullfile(user_path,'data/analysis/sample_GaSb_W');
    input_param.save_formats      = 'mat';

    mkdir(input_param.save_folder);
    

    input_param.detector_id         = 'merlin';
    scan_d = Scan(input_param);
    
    input_param.detector_id         = 'xspress3';
    scan_f = Scan(input_param);

%     input_param.detector_id         = 'pil100k';
%     scan_p = Scan(input_param);

    % Reading
    scan_d.create_path(input_param.pre_path_data);
    scan_d.read_nanomax_merlin;

    scan_f.create_path(input_param.pre_path_data);
    scan_f.read_nanomax_xspress3;
    
%     scan_p.create_path(input_param.pre_path_data);
%     scan_p.read_nanomax_pil100k;
    
        % Corrections
    scan_d.read_mask;
    scan_d.correct_mask_nanomax;
    scan_d.correct_flux;

%      scan_p.correct_flux; 
    
    % Saving
    scan_d.data = permute(scan_d.data,[2,1,3,4]); 
    scan_d.flip('lr'); 
    scan_d.save_data;

    scan_f.save_data; 
    
%     scan_p.save_data;
    
    disp('### Data import done! ###');
end
