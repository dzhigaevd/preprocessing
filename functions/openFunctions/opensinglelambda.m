function out = opensinglelambda(varargin)
% ---
% --- OPENSINGLELAMBDA Open LAMBDA hdf5 (*.h5) or NEXUS (*.nxs) data files
% ---
% --- USAGE:
% --- OUT = OPENSINGLELAMBDA(FILENAME)        --> direct call to a file
% --- OUT = OPENSINGLELAMBDA(FILENAME,INDEX)  --> call to open a specific 
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
% --- $Revision: 1.0 $Date: 2015/05/21 $ function to open LAMBDA hdf5 files
% --- $Revision: 1.1 $Date: 2015/06/18 $ function to open LAMBDA nxs files
% ---

% =========================================================================
% --- At a certain point, Lambda images will be stored distortion corrected
% =========================================================================
isDistCorrected      = 0                                                   ; % Switch to indicate if the data was stored distortion corrected
DistortionCorrection = datenum(datestr(datevec('2015-07-01','yyyy-mm-dd'))); % MS 20150619 (to be finalized)
isHeaderCorrected    = 0                                                   ; % Switch to indicate if the NEXUS header is in new format
HeaderCorrection     = datenum(datestr(datevec('2016-04-01','yyyy-mm-dd'))); % MS 20160407

% =========================================================================
% --- The Lambda 02 detector has 'moving hot pixels' in a small region
% =========================================================================
hasMovingHotPixels = 0                                                     ; % Switch to indicate if the data contains 'moving hot pixels'

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
[~,~,ext] = fileparts(fullhfile)                                           ;
fclose(fid)                                                                ;

% =========================================================================
% --- Determine some standard parameters from the h5-file
% --- Try to read out HDF_5 header file...
% =========================================================================
% --- Retreive Start Time & End Time
MyStartTime  = h5readatt(fullhfile, '/', 'file_time' )                     ;
MyEndTime    = h5readatt(fullhfile, '/', 'file_update_time' )              ;

% --- Calculate time
if str2double(MyStartTime{1}(12)) > str2double(MyEndTime{1}(12))
    Count_hours = str2double(MyEndTime{1}(12:13)) - str2double(MyStartTime{1}(12:13)) + 24;
else
    Count_hours = str2double(MyEndTime{1}(12:13)) - str2double(MyStartTime{1}(12:13));
end
if str2double(MyStartTime{1}(15)) > str2double(MyEndTime{1}(15))
    Count_minutes = str2double(MyEndTime{1}(15:16)) - str2double(MyStartTime{1}(15:16)) + 60;
    Count_hours   = Count_hours - 1;
else
    Count_minutes = str2double(MyEndTime{1}(15:16)) - str2double(MyStartTime{1}(15:16));
end
if str2double(MyStartTime{1}(18)) > str2double(MyEndTime{1}(18))
    Count_seconds = str2double(MyEndTime{1}(18:19)) - str2double(MyStartTime{1}(18:19)) + 60;
    Count_minutes = Count_minutes - 1;
else
    Count_seconds = str2double(MyEndTime{1}(18:19)) - str2double(MyStartTime{1}(18:19));
end
if str2double(MyStartTime{1}(21)) > str2double(MyEndTime{1}(21))
    Count_mseconds = str2double(MyEndTime{1}(21:23)) - str2double(MyStartTime{1}(21:23)) + 1000;
    Count_seconds  = Count_seconds - 1;
else
    Count_mseconds = str2double(MyEndTime{1}(21:23)) - str2double(MyStartTime{1}(21:23));
end
Count_complete = (Count_hours.*3600) + (Count_minutes.*60) + Count_seconds + (Count_mseconds./1000);
clear Count_hours Count_minutes Count_seconds Count_mseconds

% =========================================================================
% --- Compare MyStartTime to the Distortion Correction Date
% =========================================================================
% --- Convert the 'file_time' variable to a date number
CompareDateStr = MyStartTime{1}                                            ;
CompareDateStr = CompareDateStr(1:10)                                      ;
CompareDate    = datenum(datestr(datevec(CompareDateStr,'yyyy-mm-dd')))    ;
% --- If needed, enable switch to indicate if data is distortion corrected
if ( CompareDate >  DistortionCorrection )
    isDistCorrected = 1                                                    ;   
end
if ( CompareDate >  HeaderCorrection )
    isHeaderCorrected = 1                                                  ;   
end
% disp(isDistCorrected)
clear MyStartTime MyEndTime CompareDateStr CompareDate                     ;
clear DistortionCorrection HeaderCorrection                                ;

% =========================================================================
% --- Open HDF5 file - Retrieve some meta data
% =========================================================================
fid = H5F.open(fullhfile, 'H5F_ACC_RDONLY', 'H5P_DEFAULT')                 ;

% --- Get metadata starting from the root group
parent_group_id    = H5G.open(fid, '/', 'H5P_DEFAULT')                     ;
parent_group_name  = H5I.get_name(parent_group_id)                         ;
child1_group_name  = '/'                                                   ;
child1_full_name   = child1_group_name;
child1_group_id    = H5G.open(parent_group_id, child1_group_name, 'H5P_DEFAULT') ;
child1_group_info  = H5G.get_info(child1_group_id)                         ;
child1_group_info  = child1_group_info.nlinks                              ;

% --- Only 1 subfolder to descend into... Link 0
child2_group_name  = H5L.get_name_by_idx(parent_group_id, child1_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', 0, 'H5P_DEFAULT') ;
child2_full_name   = [child1_full_name, child2_group_name]                 ;
child2_group_id    = H5G.open(child1_group_id, child2_group_name, 'H5P_DEFAULT') ;    
child2_group_info  = H5G.get_info(child2_group_id)                         ;
child2_group_info  = child2_group_info.nlinks                              ;

% --- Only 1 subfolder to descend into... Link 0
child3_group_name  = H5L.get_name_by_idx(child1_group_id, child2_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', 0, 'H5P_DEFAULT') ;
child3_full_name   = [child2_full_name, '/' child3_group_name]             ;
child3_group_id    = H5G.open(child2_group_id, child3_group_name, 'H5P_DEFAULT') ;
child3_group_info  = H5G.get_info(child3_group_id)                         ;
child3_group_info  = child3_group_info.nlinks                              ;

% --- Only 1 subfolder to descend into... Link 0
child4_group_name  = H5L.get_name_by_idx(child2_group_id, child3_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', 0, 'H5P_DEFAULT') ;
child4_full_name   = [child3_full_name, '/', child4_group_name]            ;
child4_group_id    = H5G.open(child3_group_id, child4_group_name, 'H5P_DEFAULT') ;
child4_group_info  = H5G.get_info(child4_group_id)                         ;
child4_group_info  = child4_group_info.nlinks                              ;

% --- Get first group name
% --- For newer Tango versions additional information are stored here..
child5_group_name = H5L.get_name_by_idx(child3_group_id, child4_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', 0, 'H5P_DEFAULT') ;

% --- Check object-type
obj_id   =  H5O.open(child4_group_id, child5_group_name, 'H5P_DEFAULT')    ;
obj_info =  H5O.get_info(obj_id)                                           ;
H5O.close(obj_id)                                                          ;

if obj_info.type == 0 % A real link is existing !
    FirstDatasets = cell(child4_group_info,1)                              ;
    FirstDatasets{1} = 'Dummy'                                             ;
    for mm = 2:child4_group_info
        FirstDatasets{mm} = H5L.get_name_by_idx(child3_group_id, child4_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', (mm-1), 'H5P_DEFAULT') ;
    end
    child5_full_name  = [child4_full_name, '/', child5_group_name]         ;
    child5_group_id   = H5G.open(child4_group_id, child5_group_name, 'H5P_DEFAULT') ;
    child5_group_info = H5G.get_info(child5_group_id)                      ;
    child5_group_info = child5_group_info.nlinks                           ;
    
    % --- Load second set of datasets
    SecondDatasets = cell(child5_group_info,1)                             ;
    for mm = 1:child5_group_info
        SecondDatasets{mm} = H5L.get_name_by_idx(child4_group_id, child5_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', (mm-1), 'H5P_DEFAULT') ;
    end
    H5G.close(child5_group_id)                                             ;
else
    FirstDatasets = cell(child4_group_info,1)                              ;
    for mm = 1:child4_group_info
        FirstDatasets{mm} = H5L.get_name_by_idx(child3_group_id, child4_group_name, 'H5_INDEX_NAME', 'H5_ITER_INC', (mm-1), 'H5P_DEFAULT') ;
    end
end

% --- Get size of dataset
obj_id   = H5O.open(child4_group_id, 'data', 'H5P_DEFAULT')                ;
space_id = H5D.get_space(obj_id)                                           ;
[~,dims] = H5S.get_simple_extent_dims(space_id)                            ;
Datasize = fliplr(dims)                                                    ;
H5O.close(obj_id)                                                          ;

H5G.close(child4_group_id)                                                 ;
H5G.close(child3_group_id)                                                 ;
H5G.close(child2_group_id)                                                 ;
H5G.close(child1_group_id)                                                 ;
H5G.close(parent_group_id)                                                 ;

Data_prefix1        = [child4_full_name, '/data/']                         ;
if strcmpi(ext,'.h5') == 1 || strcmpi(ext,'.nxs') == 1
    Data_prefix2 = [child4_full_name, '/count_time/']                      ;           
    try
        count_time = h5read(file, Data_prefix2) / 1000                     ;
    catch
        Data_prefix2 = [child4_full_name, '/collection/count_time/']       ;
        try
            count_time = h5read(file, Data_prefix2) / 1000                 ;
        catch
        end
        
    end
end

% --- Load in additional information
Data_prefix11 = [child4_full_name, '/delay_time/']                         ;
try
    Delay_time          = h5read(file, Data_prefix11)                      ;
catch
    try
        Data_prefix11 = [child4_full_name, '/collection/delay_time/']      ;
        Delay_time    = h5read(file, Data_prefix11)                        ;
    catch
    end
end
Data_prefix12 = [child4_full_name, '/frame_numbers/']                      ;
try
    Frame_Number        = h5read(file, Data_prefix12)                      ;
catch
    try 
        Data_prefix12 = [child4_full_name, '/collection/frame_numbers/']   ;
        Frame_Number  = h5read(file, Data_prefix12)                        ;
    catch
    end
end
Data_prefix13 = [child4_full_name, '/total_loss_frames/']                  ;
try
    Total_Loss          = h5read(file, Data_prefix13)                      ;
catch
    try
        Data_prefix13 = [child4_full_name, '/collection/total_loss_frames/'] ;
        Total_Loss    = h5read(file, Data_prefix13)                        ;
    catch
    end
end
Data_prefix14 = [child4_full_name, '/shutter_time/']                       ;
try
    Shutter_time = h5read(file, Data_prefix14)                             ;
catch
    try
        Data_prefix14 = [child4_full_name, '/collection/shutter_time/']    ;
        Shutter_time = h5read(file, Data_prefix14)                         ;
    catch
    end
end

% --- Check illumination time
if exist('count_time', 'var') == 0 && exist('Shutter_time', 'var') 
    count_time = Shutter_time / 1000                                       ;
end

% --- Check number of frames for consistency
if (exist('Frame_Number', 'var') == 1)
    if Datasize(3) ~= Frame_Number
        warning('Check number of frames !')
    end
end

% --- Check illumination time
if exist('count_time', 'var') == 1
    CompleteLength = Datasize(3) * count_time                              ;
    if CompleteLength > Count_complete
        warning('Complete counting time mismatch !')
    end
end

% --- Check for lost frames
if (exist('Total_Loss', 'var') == 1)
    if Total_Loss ~= 0
        % warning('Frames lost during aquisition !')
    end
end

% =========================================================================
% --- Adjusted due to additional introduced rows/columns
% =========================================================================
if isDistCorrected == 0
    rows    = 1556                                                         ;
    cols    = 516                                                          ;
else
    rows    = Datasize(1)                                                  ;    
    cols    = Datasize(2)                                                  ;
end

bytes   = 4                                                                ;
dlen    = bytes * cols * rows                                              ;
col_beg = 1                                                                ;
col_end = cols                                                             ;
row_beg = 1                                                                ;
row_end = rows                                                             ;

preset  = []                                                               ;
prefix  = ''                                                               ;
suffix  = ext                                                              ;
number  = []                                                               ;

% =========================================================================
% --- create variables, if certain expected header information were missing
% =========================================================================
if ( exist('count_time','var') == 0 )
    count_time = 1.0                                                       ;
    disp('No count_time variable')
end
if (  exist('number','var') == 0 )
    number = []                                                            ;
end
if ( exist('row_beg','var') == 0 )
    row_beg = 1                                                            ;
end
if ( exist('row_end','var') == 0 )
    row_end = rows                                                         ;
end
if ( exist('col_beg','var') == 0 )
    col_beg = 1                                                            ;
end
if ( exist('col_end','var') == 0 )
    col_end = cols                                                         ;
end
if ( exist('prefix','var') == 0 )
    prefix = ''                                                            ;
end
if ( exist('suffix','var') == 0 )
    suffix = '.h5'                                                         ;
end

% =========================================================================
% --- try to find the elapsed time between frames
% =========================================================================
preset  = count_time                                                       ;
elapsed = count_time                                                       ;

% =========================================================================
% --- create field names of standard edf header structure
% --- leave values empty
% =========================================================================

h5header(1,:)  = {'HeaderID',              []}                             ;
h5header(2,:)  = {'ByteOrder',             []}                             ;
h5header(3,:)  = {'DataType',              []}                             ;
h5header(4,:)  = {'Size',              dlen*2}                             ;
h5header(5,:)  = {'Dim_1',               rows}                             ;
h5header(6,:)  = {'Dim_2',               cols}                             ;
h5header(7,:)  = {'Image',                 []}                             ;
h5header(8,:)  = {'acq_frame_nb',          []}                             ;
h5header(9,:)  = {'time',                  []}                             ;
h5header(10,:) = {'time_of_day',           []}                             ;
h5header(11,:) = {'time_of_frame',         []}                             ;
h5header(12,:) = {'preset',            preset}                             ;

out.h5header = h5header                                                    ;

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
header(6,:)  = {'suffix',        suffix}                                   ;
header(7,:)  = {'monitor',       []}                                       ;
header(8,:)  = {'shutter',       []}                                       ;
header(9,:)  = {'row_beg',       row_beg}                                  ;
header(10,:) = {'row_end',       row_end}                                  ;
header(11,:) = {'col_beg',       col_beg}                                  ;
header(12,:) = {'col_end',       col_end}                                  ;
header(13,:) = {'row_bin',       []}                                       ;
header(14,:) = {'col_bin',       []}                                       ;
header(15,:) = {'rows',          rows}                                     ;
header(16,:) = {'cols',          cols}                                     ;
header(17,:) = {'bytes',         bytes}                                    ;
header(18,:) = {'kinetics',      []}                                       ;
header(19,:) = {'kinwinsize',    []}                                       ;
header(20,:) = {'elapsed',       elapsed}                                  ;
header(21,:) = {'preset',        preset}                                   ;
header(22,:) = {'topup',         []}                                       ;
header(23,:) = {'inject',        []}                                       ;
header(24,:) = {'dlen',          dlen / bytes}                             ; % assume it is uncompressed
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

% =========================================================================
% --- Open the file and reshape if necessary
% =========================================================================
imm            = h5read(file, Data_prefix1, [1 1 Index], [inf inf 1])      ;
imm            = single(imm)                                               ;

% --- If needed, reshape frames to insert rows/columns for big pixels
if isDistCorrected == 0 
    NewCols        = [imm(:,256), imm(:,256), imm(:,256),               ...
                      imm(:,257), imm(:,257), imm(:,257)]                  ;
    imm            = [imm(:,1:255), (NewCols ./ 3), imm(:,258:512)]        ;

    NewRows1       = [imm( 256,:); imm( 256,:); imm( 256,:);            ...
                      imm( 257,:); imm( 257,:); imm( 257,:)]               ;
    NewRows2       = [imm( 512,:); imm( 512,:); imm( 512,:);            ...
                      imm( 513,:); imm( 513,:); imm( 513,:)]               ;
    NewRows3       = [imm( 768,:); imm( 768,:); imm( 768,:);            ...
                      imm( 769,:); imm( 769,:); imm( 769,:)]               ;
    NewRows4       = [imm(1024,:); imm(1024,:); imm(1024,:);            ...
                      imm(1025,:); imm(1025,:); imm(1025,:)]               ;
    NewRows5       = [imm(1280,:); imm(1280,:); imm(1280,:);            ...
                      imm(1281,:); imm(1281,:); imm(1281,:)]               ;

    imm            = [imm(   1: 255,:); (NewRows1 ./ 3);                ...
                      imm( 258: 511,:); (NewRows2 ./ 3);                ...
                      imm( 514: 767,:); (NewRows3 ./ 3);                ...
                      imm( 770:1023,:); (NewRows4 ./ 3);                ...
                      imm(1026:1279,:); (NewRows5 ./ 3);                ...
                      imm(1282:1536,:)]                                    ;
end

% --- If needed, correct for 'moving hot pixels'
% --- These pixels located at 'even' rows 746:775 and coloums 1:255
if hasMovingHotPixels == 1
    imm(746:775,1:255) = lambda_correcthotpixel(imm(746:775, 1:255))       ;
end

% --- Store image to output structure
out.imm = imm                                                              ;

% ---
% EOF

