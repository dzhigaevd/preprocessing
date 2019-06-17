function out = opennxs(varargin)
% ---
% --- OPENNXS Distributor to open NEXUS (*.nxs) data files
% ---
% --- USAGE:
% --- OUT = OPENNXS(FILENAME)          --> direct call to a file
% --- OUT = OPENNXS(FILENAME,INDEX)    --> call to open a specific 
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
    info  = h5info(file,'/entry/instrument/detector/data')                 ;
    dummy = 1                                                              ;
catch
end

% if dummy == 0
%     try
%         info  = h5info(file,'/entry/instrument/detector/data')           ;
%         dummy = 2                                                        ;
%     catch
%     end
% end

switch dummy
    case 0
        out = []                                                           ;
    case 1
        out  = opensinglelambda(file,Index)                                ; % for a NEXUS file (Lambda style)
end
% ---
% EOF
