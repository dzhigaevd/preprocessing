function out = opensaclah5(varargin)
% ---
% --- OPENSACLAH5 Open SACLA hdf5 (*.h5) data files
% ---
% --- USAGE:
% --- OUT = OPENSACLAH5(FILENAME)        --> direct call to a file
% --- OUT = OPENSACLAH5(FILENAME,INDEX)  --> call to open a specific 
% --- image from a batch of single images. Here filename stands for the 
% --- first file of the batch and the corresponding filename is created by
% --- cutting the filename after the last underscore and adding the index 
% --- to the remaining namestring.
% ---
% --- Input Argument:
% --- FILE     : image file name
% --- Index    : image index
% ---
% --- Output Argument:
% --- OUT : structure containing header & image.
% ---
% --- Pawel Kwasniewski
% --- $Revision: 1.0 $Date: 2014/01/09 $ function to open SACLA hdf5 files
% ---

% =========================================================================
% --- input file preparation
% =========================================================================
file     = varargin{1}                                                     ;
Index    = varargin{2}                                                     ;

out.imm = h5read(file,strcat(Index,'/detector_data'))                      ;

% ---
% EOF