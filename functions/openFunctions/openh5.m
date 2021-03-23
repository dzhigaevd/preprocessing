function out = openh5(varargin)
% ---
% --- OPENH5 Distributor to open hdf5 (*.h5) data files
% ---
% --- USAGE:
% --- OUT = OPENH5(FILENAME)          --> direct call to a file
% --- OUT = OPENH5(FILENAME,INDEX)    --> call to open a specific 
% --- image from a batch of single images. Here filename stands for the 
% --- first file of the batch and the corresponding filename is created by
% --- cutting the filename after the last underscore and adding the index 
% --- to the remaining namestring.
% ---
% --- Input Argument:
% --- FILENAME : image file name
% --- Index    : image index
% ---
% --- Output Argument:
% --- OUT : structure containing header & image.
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2015/05/21 $ function to open h5 files
% ---

% =========================================================================
% --- input file preparation
% =========================================================================
file     = varargin{1}                                                     ;
Index    = varargin{2}                                                     ;

dummy = 0                                                                  ;
try
    info  = h5info(file,'/detector_data')                                  ; %#ok<*NASGU>
    dummy = 1                                                              ;
catch
end

if dummy == 0
    try
        info  = h5info(file,'/entry/instrument/detector/data')             ;
        dummy = 2                                                          ;
    catch
    end
end

if dummy == 0
    try
        info  = h5info(file,'/entry/data/data')                            ;
        dummy = 3                                                          ;
    catch
        try
            info  = h5info(file,'/entry/data/data_000001')                 ;
            dummy = 3                                                      ;
        catch
        end
    end
end

switch dummy
    case 0
        out = []                                                           ;
    case 1
        out  = opensaclah5(file,Index)                                     ; % for a hdf5 (SACLA style)
    case 2
        out  = opensinglelambda(file,Index)                                ; % for a hdf5 (Lambda style)
    case 3
        out  = opensingleeiger4m(file,Index)                               ; % for a hdf5 (EIGER 4M style)
end

% ---
% EOF
