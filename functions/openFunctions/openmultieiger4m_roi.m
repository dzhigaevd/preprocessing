function out = openmultieiger4m_roi(varargin)
% ---
% --- OPENMULTIEIGER4M Open EIGER 4M hdf5 (.h5) data files
% ---
% --- USAGE
% --- OUT = OPENMULTIEIGER4M(FILENAME)        -- direct call to 
% --- open all images of a HDF5 file
% --- OUT = OPENMULTIEIGER4M(FILENAME,STARTINDEX,ENDINDEX)-- call to 
% --- open a block of images. 
% --- OUT = OPENMULTIEIGER4M(FILENAME,STARTINDEX,ENDINDEX,ROI)-- call to 
% --- open an ROI of a block of images. 
% --- OUT = OPENMULTIEIGER4M(FILENAME,ROI)-- call to 
% --- open an ROI of all images. 
% ---
% --- Input Argument
% --- FILE        image file name
% --- STARTINDEX  starting image index
% --- ENDINDEX    ending image index
% --- ROI         [roiYstart, roiYend, roiXstart, roiXend]
% ---
% --- Output Argument
% --- OUT  structure containing header & image block.
% ---
% --- M. Sprung
% --- $Revision 1.0 $Date 20160912 $ function to open EIGER 4M files
% ---

out = []                                                                   ;

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file               = varargin{1}                                       ;
    StartIndex         = 1                                                 ;
    EndIndex           = inf                                               ;
    roiYstart          = 1                                                 ;
    roiXstart          = 1                                                 ;
    roiYrange          = inf                                               ;
    roiXrange          = inf                                               ;
elseif ( nargin == 2 )
    file               = varargin{1}                                       ;
    StartIndex         = 1                                                 ;
    EndIndex           = inf                                               ;
    roiYstart          = varargin{2}(1)                                    ;
    roiYend            = varargin{2}(2)                                    ;
    roiXstart          = varargin{2}(3)                                    ;
    roiXend            = varargin{2}(4)                                    ;
    roiYrange          = roiYend - roiYstart + 1                           ;
    roiXrange          = roiXend - roiXstart + 1                           ;
    clear roiXend roiYend                                                  ;
elseif ( nargin == 3 )
    file               = varargin{1}                                       ;
    StartIndex         = varargin{2}                                       ;
    EndIndex           = varargin{3}                                       ;
    if StartIndex >= EndIndex
        dummy      = EndIndex                                              ;
        EndIndex   = StartIndex                                            ;
        StartIndex = dummy                                                 ;
        clear dummy                                                        ;
    end
    roiYstart          = 1                                                 ;
    roiXstart          = 1                                                 ;
    roiYrange          = inf                                               ;
    roiXrange          = inf                                               ;
elseif ( nargin == 4 )
    file               = varargin{1}                                       ;
    StartIndex         = varargin{2}                                       ;
    EndIndex           = varargin{3}                                       ;
    if StartIndex >= EndIndex
        dummy      = EndIndex                                              ;
        EndIndex   = StartIndex                                            ;
        StartIndex = dummy                                                 ;
        clear dummy                                                        ;
    end
    roiYstart          = varargin{4}(1)                                    ;
    roiYend            = varargin{4}(2)                                    ;
    roiXstart          = varargin{4}(3)                                    ;
    roiXend            = varargin{4}(4)                                    ;
    roiYrange          = roiYend - roiYstart + 1                           ;
    roiXrange          = roiXend - roiXstart + 1                           ;
    clear roiXend roiYend                                                  ;
end

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
        masterstart = name(1:end-12)                                        ;
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
% --- if it points to a master file, than look for the first data file
% =========================================================================
if data == 0
    if numel(name) > 7
        datastart = name(1:end-7)                                          ;
    else
        disp('Error! It is not possible to recontruct the file name to a data file') ;
    end
    datafile1           = [pathstr filesep datastart '_data_'           ...
                           num2str(1,'%06i') ext]                          ; % build the name of the first file
    [data1fid,~]        = fopen(datafile1,'r')                             ;
    if ( data1fid >= 3 )                                                     % if data file 1 opens, than set the data1 switch to 1
        data1 = 1                                                          ;
    else
        disp('Warning! It is not possible to access the first data file')
    end
    try
        fclose(data1fid)                                                   ;
    catch
    end
end

% =========================================================================
% --- Retrieve and reshape the data if necessary
% --- Multiple cases are possible
% --- Case 1 'file' points to a 'data' file
% --- Case 2 'file' points to a 'master' file (here will be sub cases)
% =========================================================================
% --- Case 1 'file' points directly to a 'data' file
if data == 1
    dataInfo = h5info(file, '/entry/data/data/')                           ;
    if EndIndex == inf
        EndIndex = dataInfo.Dataspace.Size(3)                              ; % set EndIndex to the chunk size
    end
    if dataInfo.Dataspace.Size(3) >= StartIndex && dataInfo.Dataspace.Size(3) >= EndIndex
        for j = 1 : EndIndex - StartIndex + 1
            img = h5read(file,'/entry/data/data/',[roiXstart roiYstart StartIndex+j-1],[roiXrange roiYrange 1]) ;
            if j == 1
                try
                    out.imm = zeros(size(img,2),size(img,1),EndIndex-StartIndex+1,'single') ;
                catch EM
                    disp(EM)                                               ;
                    fprintf('Please load the data in smaller chunks!')     ;
                    return                                                 ;
                end
            end
            out.imm(:,:,j) = single(img.')                                 ;
            if mod(j,500) == 0
                disp(j)
            end
        end
    else
        fprintf('Error! Data file contains only %i frames!',dataInfo.Dataspace.Size(3)) ; %DSPS
        return
    end
% --- Case 2 'file' points to a 'master' file
elseif master == 1
    if data1 == 1                                                            % (at least the first data file exists)
        dataInfo = h5info(datafile1, '/entry/data/data/')                  ;
        
        % --- the first data file contains all the 'wanted' images
        if dataInfo.Dataspace.Size(3) >= StartIndex                     ...
        && dataInfo.Dataspace.Size(3) >= EndIndex                          
            for j = 1 : EndIndex - StartIndex + 1
                img = h5read(datafile1,'/entry/data/data/',[roiXstart roiYstart StartIndex+j-1],[roiXrange roiYrange 1]) ;
                if j == 1
                    out.imm = zeros(size(img,2),size(img,1),EndIndex-StartIndex+1,'single')  ; 
                end
                out.imm(:,:,j) = single(img.')                             ;
                if mod(j,500) == 0
                    disp(j)
                end
            end
        % --- the 'wanted' images are possibly in different 'data' file(s)
        else
            
            % --- Calculate total number of image pointers in master file
            if EndIndex == inf
                dummybreak  = 0                                            ;
                lastdatanum = 0                                            ;
                EndIndex    = 0                                            ; % reset 'EndIndex'
                while lastdatanum ~= -1 && dummybreak < 10000
                    dummybreak      = dummybreak + 1                       ;
                    lastdatanum     = lastdatanum + 1                      ;
                    lastdatafile    = [pathstr filesep datastart '_data_'  ...
                                       num2str(lastdatanum,'%06i') ext]    ; % build the file name to the 'next' data set
                    [lastdatafid,~] = fopen(lastdatafile,'r')              ;
                    if ( lastdatafid == -1 )                                 % if Matlab can't open the 'lastdatafile' it returns '-1'
                        lastdatanum = -1                                   ; % create the 'break' condition
                    elseif (lastdatafid >= 3)
                        lastdataInfo = h5info(lastdatafile,'/entry/data/data/')  ;
                        EndIndex     = EndIndex + lastdataInfo.Dataspace.Size(3) ; 
                    end
                    try
                        fclose(lastdatafid)                                ;
                    catch
                    end
                end
            end
            
            goal1    = ceil(StartIndex / dataInfo.Dataspace.Size(3))       ; % data file number of StartIndex
            goal2    = ceil(EndIndex   / dataInfo.Dataspace.Size(3))       ; % data file number of EndIndex
            
            % --- the 'wanted' images resides in a single data file
            if goal1 == goal2
                newStartIndex = StartIndex - (goal1-1) * dataInfo.Dataspace.Size(3) ;
                newEndIndex   = EndIndex   - (goal2-1) * dataInfo.Dataspace.Size(3) ;
                
                datafilegoal  = [pathstr filesep datastart '_data_'     ...
                                 num2str(goal1,'%06i') ext]                ; % build the file name to the 'goal' data set
                [datagoalfid,~] = fopen(datafilegoal,'r')                  ;
                if ( datagoalfid < 3 )
                    disp('Error: Goal data file can not be opened!')       ;
                end
                try
                    fclose(datagoalfid)                                    ;
                catch
                end
                dataInfoGoal = h5info(datafilegoal, '/entry/data/data/')   ;
                
                if dataInfoGoal.Dataspace.Size(3) >= newStartIndex      ...
                && dataInfoGoal.Dataspace.Size(3) >= newEndIndex            % the goal data file contains the 'wanted' image
                    for j = 1 : newEndIndex - newStartIndex + 1
                        img = h5read(datafilegoal,'/entry/data/data/',[roiXstart roiYstart j],[roiXrange roiYrange 1]) ;
                        if j == 1
                            try
                                out.imm = zeros(size(img,2),size(img,1),newEndIndex-newStartIndex+1,'single')  ;
                            catch EM
                                disp(EM)                                   ;
                                return                                     ;
                            end
                        end
                        out.imm(:,:,j) = single(img.')                     ;
                        if mod(j,500) == 0
                            disp(j)
                        end
                    end
                else                                                        % the 'wanted' images are in a multiple 'data' file
                    fprintf('Error: Goal data file contains only %i frames!',dataInfoGoal.Dataspace.Size(3)) ;
                    return
                end
                
            % --- the 'wanted' images resides in a neighboring data files
            elseif (goal2 - goal1) == 1
                
                newStartIndex = StartIndex - (goal1-1) * dataInfo.Dataspace.Size(3) ;
                newEndIndex   = EndIndex   - (goal2-1) * dataInfo.Dataspace.Size(3) ;
                
                datafilegoal1  = [pathstr filesep datastart '_data_'     ...
                                 num2str(goal1,'%06i') ext]                ; % build the file name to the 'goal1' data set
                [datagoal1fid,~] = fopen(datafilegoal1,'r')                ;
                if ( datagoal1fid < 3 )
                    disp('Error: Goal1 data file can not be opened!')      ;
                end
                try
                    fclose(datagoal1fid)                                   ;
                catch
                end
                dataInfoGoal1 = h5info(datafilegoal1, '/entry/data/data/') ;
                
                datafilegoal2  = [pathstr filesep datastart '_data_'     ...
                                 num2str(goal2,'%06i') ext]                ; % build the file name to the 'goal2' data set
                [datagoal2fid,~] = fopen(datafilegoal2,'r')                ;
                if ( datagoal2fid < 3 )
                    disp('Error: Goal2 data file can not be opened!')      ;
                end
                try
                    fclose(datagoal2fid)                                   ;
                catch
                end
                dataInfoGoal2 = h5info(datafilegoal2, '/entry/data/data/') ;
                
                if dataInfoGoal1.Dataspace.Size(3) >= newStartIndex     ...
                && dataInfoGoal2.Dataspace.Size(3) >= newEndIndex            % the goal data sets contain the 'wanted' images
                    for k = newStartIndex : dataInfoGoal1.Dataspace.Size(3)
                        img = h5read(datafilegoal1,'/entry/data/data/',[roiXstart roiYstart k],[roiXrange roiYrange 1]) ;
                        if k == newStartIndex
                            try
                                out.imm = zeros(size(img,2),size(img,1),EndIndex-StartIndex+1,'single')  ;
                            catch EM
                                disp(EM)                                   ;
                                return                                     ;
                            end
                        end
                        j = k - newStartIndex + 1                          ;
                        out.imm(:,:,j) = single(img.')                     ;
                    end
                    for k = 1 : newEndIndex
                        img = h5read(datafilegoal2,'/entry/data/data/',[roiXstart roiYstart k],[roiXrange roiYrange 1]) ;
                        j = j + 1                                          ;
                        out.imm(:,:,j) = single(img.')                     ;
                    end
                else
                    fprintf('Error: Images could not be retrieved from data sets!') ;
                    return
                end
                                
            % --- images from a series of more than 2 datasets
            else     
                
                newStartIndex = StartIndex - (goal1-1) * dataInfo.Dataspace.Size(3) ;
                newEndIndex   = EndIndex   - (goal2-1) * dataInfo.Dataspace.Size(3) ;
                
                datafilegoal1  = [pathstr filesep datastart '_data_'     ...
                                 num2str(goal1,'%06i') ext]                ; % build the file name to the 'goal1' data set
                [datagoal1fid,~] = fopen(datafilegoal1,'r')                ;
                if ( datagoal1fid < 3 )
                    disp('Error: Goal1 data file can not be opened!')      ;
                end
                try
                    fclose(datagoal1fid)                                   ;
                catch
                end
                dataInfoGoal1 = h5info(datafilegoal1, '/entry/data/data/') ;
                
                datafilegoal2  = [pathstr filesep datastart '_data_'     ...
                                  num2str(goal2,'%06i') ext]               ; % build the file name to the 'goal2' data set
                [datagoal2fid,~] = fopen(datafilegoal2,'r')                ;
                if ( datagoal2fid < 3 )
                    disp('Error: Goal2 data file can not be opened!')      ;
                end
                try
                    fclose(datagoal2fid)                                   ;
                catch
                end
                dataInfoGoal2 = h5info(datafilegoal2, '/entry/data/data/') ;
                
                datafilegoaln = cell(1, goal2-goal1-1)                     ;
                dataInfoGoaln = zeros(goal2-goal1 -1, 3)                   ;
                for n = 1 : goal2 - goal1 - 1
                    datafilegoaln{n} = [pathstr filesep datastart '_data_' ...
                                        num2str(goal1 + n,'%06i') ext]     ; % build the file name to data set 'n'
                    [datagoalnfid,~] = fopen(datafilegoal2,'r')            ;
                    if ( datagoalnfid < 3 )
                        disp('Error: Goal n data file can not be opened!') ;
                    end
                    try
                        fclose(datagoalnfid)                               ;
                    catch
                    end
                    clear datagoalnfid
                    temp = h5info(datafilegoaln{n}, '/entry/data/data/')   ;
                    dataInfoGoaln(n,:) = temp.Dataspace.Size               ;
                end
                
                if dataInfoGoal1.Dataspace.Size(3) >= newStartIndex     ...
                && dataInfoGoal2.Dataspace.Size(3) >= newEndIndex            % the goal data sets contain the 'wanted' images
                    for k = newStartIndex : dataInfoGoal1.Dataspace.Size(3)
                        img = h5read(datafilegoal1,'/entry/data/data/',[roiXstart roiYstart k],[roiXrange roiYrange 1]) ;
                        if k == newStartIndex
                            try
                                out.imm = zeros(size(img,2),size(img,1),EndIndex-StartIndex+1,'single')  ;
                            catch EM
                                disp(EM)                                   ;
                                return                                     ;
                            end
                        end
                        j = k - newStartIndex + 1                          ;
                        out.imm(:,:,j) = single(img.')                     ;
                        if mod(k,500) == 0
                            disp(j)
                        end
                    end
                    for n = 1 : goal2 - goal1 - 1
                        for k = 1 : dataInfoGoaln(n,3)
                            img = h5read(datafilegoaln{n},'/entry/data/data/',[roiXstart roiYstart k],[roiXrange roiYrange 1]) ;
                            j = j + 1                                      ;
                            out.imm(:,:,j) = single(img.')                 ;
                            if mod(k,500) == 0
                                disp(j)
                            end
                        end
                    end
                    for k = 1 : newEndIndex
                        img = h5read(datafilegoal2,'/entry/data/data/',[roiXstart roiYstart k],[roiXrange roiYrange 1]) ;
                        j = j + 1                                          ;
                        out.imm(:,:,j) = single(img.')                     ;
                        if mod(k,500) == 0
                            disp(j)
                        end
                    end
                else
                    fprintf('Error: Images could not be retrieved from data sets!') ;
                    return
                end
                
            end
        end
    else
        disp('Error: It is not possible to retrieve the number of images per data file')
        return
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
        elapsed = h5read(masterfile,'/entry/instrument/detectorframe_time');
    catch
    end
    try
        preset  = h5read(masterfile,'/entry/instrument/detectorcount_time');
    catch
    end
    try
        frequency = 1 / h5read(masterfile,'/entry/instrument/detector/detectorSpecific/frame_count_time') ;
    catch
    end
end

% =========================================================================
% --- Mask bad (hot & cold) pixel as well as stripes between modules
% =========================================================================
out.imm(out.imm > 2 * 2^bitdepth * frequency * preset) = 0                 ; % maximum possible count rate for this exposure time
% out.imm(out.imm > 2^32 -  2^16 - 1) = 0                                    ;

% =========================================================================
% --- create field names of standard imm header structure
% --- leave values empty
% =========================================================================
header       = cell(53,2)                                                  ; % initialize '.header' structure
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
% --- new fields (added 102006)
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
