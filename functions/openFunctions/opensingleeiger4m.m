function out = opensingleeiger4m(varargin)
% ---
% --- OPENSINGLEEIGER4M Open EIGER 4M hdf5 (*.h5)data files
% ---
% --- USAGE:
% --- OUT = OPENSINGLEEIGER4M(FILENAME)        --> direct call to a file
% --- OUT = OPENSINGLEEIGER4M(FILENAME,INDEX)  --> call to open a specific 
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
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2016/09/12 $ function to open EIGER 4M files
% ---

out = []                                                                   ;

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to the first frame of a nexus file
    Index    = 1                                                           ;
elseif ( nargin == 2 )
    file               = varargin{1}                                       ;
    Index              = varargin{2}                                       ;
end

% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
fullhfile = fopen(fid)                                                     ;  % get full pathname
[pathstr,name,ext] = fileparts(fullhfile)                                  ;
fclose(fid)                                                                ;

% =========================================================================
% --- check if 'file' points to a master or a data file
% =========================================================================
master = 0                                                                 ; % switch for master file
data   = 0                                                                 ; % switch for data file
data1  = 0                                                                 ; % switch for first data file
if numel(name) > 7 && strcmpi(name(end-6:end),'_master') == 1
    master     = 1                                                         ;
    masterfile = file                                                      ;
elseif numel(name) > 12 && strcmpi(name(end-11:end-6),'_data_') == 1
    data = 1                                                               ;
else
    disp('Error! Input File not recognized as master or data file')        ;
    return
end

% =========================================================================
% --- if it points to a data file, than look for the master file
% =========================================================================
if master == 0
    if numel(name) > 12
        masterstart = name(1:end-12)                                       ;
    end
    masterfile          = [pathstr filesep masterstart '_master' ext]      ; % build the name of the master file
    [masterfid,~]       = fopen(masterfile,'r')                            ;
    if ( masterfid >= 3 )                                                    % if master file opens, than set the master switch to 1
        master = 1                                                         ;
    end
    try
        fclose(masterfid)                                                  ;
    catch
    end
end

% =========================================================================
% --- if it points to a master file, than look for data file(s)
% =========================================================================
if data == 0
    if numel(name) > 7
        datastart = name(1:end-7)                                          ;
    else
        disp('Error! It is not possible to recontruct the file name of a data file') ;
    end
    datafile1           = [pathstr filesep datastart '_data_'           ...
                           num2str(1,'%06i') ext]                          ; % build the name of the first file
    [data1fid,~]        = fopen(datafile1,'r')                             ;
    if ( data1fid >= 3 )                                                     % if data file 1 opens, than set the master switch to 1
        data1 = 1                                                          ;
    end
    try
        fclose(data1fid)                                                   ;
    catch
    end
end

% =========================================================================
% --- Retrieve and reshape the data if necessary
% =========================================================================
if data == 1                                                                 % file points directly to a 'data' file
    dataInfo = h5info(file, '/entry/data/data/')                           ;
    if dataInfo.Dataspace.Size(3) >= Index
        imm      = h5read(file, '/entry/data/data/', [1 1 Index], [inf inf 1]) ;
        imm      = single(imm.')                                           ;
    else
        disp(sprintf('Error! Data file contains only %i frames!',dataInfo.Dataspace.Size(3))) ; %#ok<*DSPS>
        return
    end
elseif master == 1                                                           % file points to a 'master' file (!!! data=0 !!!)
    if data1 == 1                                                            % (at least the first data file exists
        dataInfo = h5info(datafile1, '/entry/data/data/')                  ;
        if dataInfo.Dataspace.Size(3) >= Index                               % the first data file contains the 'wanted' image
            imm  = h5read(datafile1, '/entry/data/data/', [1 1 Index], [inf inf 1]) ;
            imm  = single(imm.')                                           ;
        else                                                                 % the 'wanted' image is in a different 'data' file
            goal     = ceil(Index / dataInfo.Dataspace.Size(3))            ;
            newIndex = Index - (goal-1) * dataInfo.Dataspace.Size(3)       ;
            datafilegoal    = [pathstr filesep datastart '_data_'       ...
                               num2str(goal,'%06i') ext]                   ; % build the file name to the 'goal' data set
            [datagoalfid,~] = fopen(datafilegoal,'r')                      ;
            if ( datagoalfid < 3 )
                disp('Error: Goal data file can not be opened!')           ;
            end
            try
                fclose(datagoalfid)                                        ;
            catch
            end
            dataInfoGoal = h5info(datafilegoal, '/entry/data/data/')       ;
            if dataInfoGoal.Dataspace.Size(3) >= newIndex                    % the first data file contains the 'wanted' image
                imm  = h5read(datafilegoal, '/entry/data/data/', [1 1 newIndex], [inf inf 1]) ;
                imm  = single(imm.')                                       ;
            else                                                             % the 'wanted' image is in a different 'data' file
                disp(sprintf('Error: Goal data file contains only %i frames!',dataInfoGoal.Dataspace.Size(3))) ;
                return
                % case 1: goal2 - goal1 = 1
                % case 2: goal2 - goal1 > 1
            end
        end
    end
end

% =========================================================================
% --- Define 'dummy' values for missing variables
% =========================================================================
prefix    = ''                                                             ;
number    = []                                                             ;
elapsed   = 1.0                                                            ;
preset    = 1.0                                                            ;
frequency = 760.0                                                          ;
bitdepth  = 12                                                             ;
if master == 1
    try
        elapsed = h5read(masterfile,'/entry/instrument/detector/frame_time/') ;
    catch
    end
    try
        preset  = h5read(masterfile,'/entry/instrument/detector/count_time/') ;
    catch
    end
    try
        frequency = 1 / h5read(masterfile,'/entry/instrument/detector/detectorSpecific/frame_count_time') ;
    catch
    end
end

% =========================================================================
% --- Mask bad (hot & cold) pixel as well as the stripes between modules
% =========================================================================
imm(imm > 2 * 2^bitdepth * frequency * preset) = 0                         ;
% imm(imm > 2^32 -  2^16 - 1) = 0                                            ;

% =========================================================================
% --- Store image to output structure
% =========================================================================
out.imm = imm                                                              ;

% =========================================================================
% --- create field names of standard imm header structure
% --- leave values empty
% =========================================================================
header = cell(53,2)                                                        ; % initialize '.header' structure
header(1,:)  = {'mode',          []}                                       ;
header(2,:)  = {'compression',   []}                                       ; 
header(3,:)  = {'date',          []}                                       ;
header(4,:)  = {'prefix',        prefix}                                   ;
header(5,:)  = {'number',        number}                                   ;
header(6,:)  = {'suffix',        ext}                                      ;
header(7,:)  = {'monitor',       []}                                       ;
header(8,:)  = {'shutter',       []}                                       ;
header(9,:)  = {'row_beg',       1}                                        ;
header(10,:) = {'row_end',       size(out.imm,1)}                          ;
header(11,:) = {'col_beg',       1}                                        ;
header(12,:) = {'col_end',       size(out.imm,2)}                          ;
header(13,:) = {'row_bin',       []}                                       ;
header(14,:) = {'col_bin',       []}                                       ;
header(15,:) = {'rows',          size(out.imm,1)}                          ;
header(16,:) = {'cols',          size(out.imm,2)}                          ;
header(17,:) = {'bytes',         4}                                        ;
header(18,:) = {'kinetics',      []}                                       ;
header(19,:) = {'kinwinsize',    []}                                       ;
header(20,:) = {'elapsed',       elapsed}                                  ;
header(21,:) = {'preset',        preset}                                   ;
header(22,:) = {'topup',         []}                                       ;
header(23,:) = {'inject',        []}                                       ;
header(24,:) = {'dlen',          4*numel(out.imm)}                         ; % assume it is uncompressed
header(25,:) = {'roi_number',    []}                                       ;
% --- modified position as of 20060306
header(26,:) = {'buffer_number', []}                                       ;
header(27,:) = {'systick',       0}                                        ;
% --- shifted header positions as of 20060306
header(28,:) = {'pv1',           []}                                       ;
header(29,:) = {'pv1VAL',        []}                                       ;
header(30,:) = {'pv2',           []}                                       ;
header(31,:) = {'pv2VAL',        []}                                       ;
header(32,:) = {'pv3',           []}                                       ;
header(33,:) = {'pv3VAL',        []}                                       ;
header(34,:) = {'pv4',           []}                                       ;
header(35,:) = {'pv4VAL',        []}                                       ;
header(36,:) = {'pv5',           []}                                       ;
header(37,:) = {'pv5VAL',        []}                                       ;
header(38,:) = {'pv6',           []}                                       ;
header(39,:) = {'pv6VAL',        []}                                       ;
header(40,:) = {'pv7',           []}                                       ;
header(41,:) = {'pv7VAL',        []}                                       ;
header(42,:) = {'pv8',           []}                                       ;
header(43,:) = {'pv8VAL',        []}                                       ;
header(44,:) = {'pv9',           []}                                       ;
header(45,:) = {'pv9VAL',        []}                                       ;
header(46,:) = {'pv10',          []}                                       ;
header(47,:) = {'pv10VAL',       []}                                       ;
% --- new fields (added 10/2006)
header(48,:) = {'imageserver',   []}                                       ;
header(49,:) = {'CPUspeed',      []}                                       ;
header(50,:) = {'immversion',    []}                                       ;
header(51,:) = {'corecotick',    0}                                        ;
header(52,:) = {'cameratype',    []}                                       ;
header(53,:) = {'threshhold',    []}                                       ;


% =========================================================================
% --- store header to output structure
% =========================================================================
out.header = header                                                        ;
clear header                                                               ;

% ---
% EOF
