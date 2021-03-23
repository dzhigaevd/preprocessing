function out = openh5attribute(varargin)

out = [];       

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
file               = varargin{1}                                           ;
attribute          = varargin{2}                                           ;

% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
fullhfile          = fopen(fid)                                            ;  % get full pathname
[pathstr,name,ext] = fileparts(fullhfile)                                  ;
fclose(fid)                                                                ;

% =========================================================================
% --- check if 'file' points to a master or a data file
% =========================================================================
master = 0                                                                 ; % switch for master file
data   = 0                                                                 ; % switch for data file
data1  = 0                                                                 ; % switch for first data file
if strcmp(ext,'.h5')
    master     = 1                                                         ;
    masterfile = file                                                      ;
    disp('Using master file!');
elseif strcmp(ext,'.hdf5')
    data = 1                                                               ;
    datafile = file;
    disp('Using data file!');
else
    disp('Error! Input File not recognized as master or data file')        ;
    return
end

% dataInfo = h5info(file, attribute);
% dataFullSize  = dataInfo.Dataspace.Size;

out = h5read(file,attribute);

%This is the case for Ni6602_buff, not checked if it is the same for more attributes
out = transpose(out);
out = out(:,any(out));