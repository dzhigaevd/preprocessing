function  out  = createPath_aps( input_param )
    % For 34 idc output is a list of tiff files
    % =========================================================================
    % --- Assemble 'master' file names
    % =========================================================================
    if input_param.maindirswitch == 1
        main_dir = fullfile('gpfs','current','raw'); % E.g. on Beamline Linux PCs
    
    elseif input_param.maindirswitch == 2
        if strcmp(input_param.online,'true') || input_param.online
            main_dir = fullfile('T:','current','raw'); % E.g. on Beamline Windows PCs        
        else
            main_dir = fullfile('userPath',input_param.beamtime_id, 'raw');
        end
        
    elseif input_param.maindirswitch == 3
           main_dir = fullfile('T:','p10','2018','data', input_param.beamtime_id, 'raw'); % E.g. on Office Windows PCs            ; % E.g. on  Windows PCs        
    
    elseif input_param.maindirswitch == 4 % User defined
        try input_param.pre_path % 34IDC option
           main_dir =  fullfile(input_param.pre_path,input_param.beamtime_id,['AD' input_param.beamtime_id '_' input_param.sample_name]) ;
        catch
           main_dir = fullfile('/asap3','petra3','gpfs','p10','2018','data', input_param.beamtime_id, 'raw'); % E.g. on Office Linux PCs            ; % E.g. on  Windows PCs        
        end
    end
    
    if strcmpi(input_param.sample_name(end),'_') ~= 1
        input_param.sample_name = [input_param.sample_name '_'];
    end

    master_folder = fullfile(main_dir, [input_param.beamtime_id '_' input_param.sample_name 'S' sprintf('%04i',input_param.scan_number)]);
    t          = dir(fullfile(master_folder,'*.tif'));
    
    for ii = 1:length(t)
        out(ii).file_name        = fullfile(master_folder,t(ii).name);
    end
%     fiofile = fullfile(main_dir, [input_param.sample_name sprintf('%05i',input_param.scan_number)],  ...
%                   [input_param.sample_name sprintf('%05i',input_param.scan_number) '.fio']);
end