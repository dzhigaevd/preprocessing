function [ masterfile, fiofile ] = createPath( inputParam )
    % =========================================================================
    % --- Assemble 'master' file names
    % =========================================================================
    if inputParam.maindirswitch == 1
        maindir = fullfile('gpfs','current','raw'); % E.g. on Beamline Linux PCs
    
    elseif inputParam.maindirswitch == 2
        if strcmp(inputParam.online,'true') || inputParam.online
            maindir = fullfile('T:','current','raw'); % E.g. on Beamline Windows PCs        
        else
            maindir = fullfile('userPath',inputParam.beamtimeId, 'raw');
        end
        
    elseif inputParam.maindirswitch == 3
           maindir = fullfile('T:','p10','2018','data', inputParam.beamtimeId, 'raw'); % E.g. on Office Windows PCs            ; % E.g. on  Windows PCs        
    
    elseif inputParam.maindirswitch == 4 % User defined
           maindir = fullfile('/asap3','petra3','gpfs','p10','2018','data', inputParam.beamtimeId, 'raw'); % E.g. on Office Linux PCs            ; % E.g. on  Windows PCs        
    end
    if strcmpi(inputParam.sampleName(end),'_') ~= 1
        inputParam.sampleName = [inputParam.sampleName '_']                                            ;
    end

    masterfile = fullfile(maindir, [inputParam.sampleName sprintf('%05i',inputParam.scanNumber)],  ...
                  inputParam.detector, [inputParam.sampleName sprintf('%05i',inputParam.scanNumber) '_master.h5']);

    fiofile = fullfile(maindir, [inputParam.sampleName sprintf('%05i',inputParam.scanNumber)],  ...
                  [inputParam.sampleName sprintf('%05i',inputParam.scanNumber) '.fio']);
end