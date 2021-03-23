function masterfile = createPath_single( inputParam )
    % =========================================================================
    % --- Assemble 'master' file names
    % =========================================================================
    if inputParam.maindirswitch == 1
        maindir = fullfile('gpfs','current','raw')                             ; % E.g. on Beamline Linux PCs
    elseif inputParam.maindirswitch == 2
        if strcmp(inputParam.online,'true')
            maindir = fullfile('T:','current','raw')                           ; % E.g. on Beamline Windows PCs        
        else
            maindir = fullfile('Y:','P10','2018','data',inputParam.beamtimeId, 'raw')           ;
    %         maindir = fullfile('U:','2017','data',beamtimeId, 'raw')           ;
        end
    elseif inputParam.maindirswitch == 3
           maindir = ['T:\2017\data\' inputParam.beamtimeId '\raw\']                      ; % E.g. on Office Windows PCs            ; % E.g. on  Windows PCs        
    end

%     if strcmpi(inputParam.sampleName(end),'_') ~= 1
%         inputParam.sampleName = [inputParam.sampleName '_']                                            ;
%     end

    masterfile = fullfile(maindir, [inputParam.sampleName],  ...
                  inputParam.detector, [inputParam.sampleName '_take_' sprintf('%05i',inputParam.masterNumber) '_master.h5']);
end

