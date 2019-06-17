function out = opensingleedf(varargin)
% ---
% --- OPENSINGLEEDF Open *.edf single files
% ---
% --- USAGE : 
% --- OUT = OPENSINGLEEDF(FILENAME)        --> direct call to a file
% --- OUT = OPENSINGLEEDF(FILENAME,INDEX)  --> call to open a specific 
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
% --- OUT : structure containing edfheader, header & image. Usually image 
% ---       is an intensity image with 0 for black, saturation point for 
% ---       white (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2007/11/26 $ function to open either compressed
% ---                                    or uncompresses *.edf single files
% --- $Revision: 2.0 $Date: 2013/02/15 $ bug fix for header length
% ---
% --- $Revision: 3.0 $Date: 2014/11/05 $ bug fix for data type
% ---

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to a single EDF file
elseif ( nargin == 2 )
    firstbatchimage    = varargin{1}                                       ;
    Index              = varargin{2}                                       ;
    [pathstr,name,ext] = fileparts(firstbatchimage)                        ;
    namestring         = name(1:end-4)                                     ; % remove the number part behind the last underscore
    namenumber         = sprintf('%04i',Index)                             ; % create new 4 digit number part of the name string
    file               = fullfile(pathstr,[namestring,namenumber,ext])     ; % create fullfile name of wanted image
end


% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end


% =========================================================================
% --- load lines from edf header
% =========================================================================
noheaderlines = 0                                                          ;
while 1
    headerline = fgetl(fid)                                                ; % get the next line of the header    
    check_startofheader = strfind(headerline,'{')                          ; % search for the start of the header
    if isempty(check_startofheader)                                          % if this is not the first line of a header continue
        check_endofheader   = strfind(headerline,'}')                      ; % search for the end of the header
        if ~isempty(check_endofheader)                                       % if the end of the header is reached break this loop
            break
        else
            noheaderlines = noheaderlines + 1                              ;
            equal_pos     = find(headerline == '=')                        ; % each line needs to contain an equal sign
            semicolon_pos = find(headerline == ';')                        ; % each line needs to be terminated by a semicolon
            right_str = headerline(equal_pos+1:semicolon_pos-1)            ;
            left_str  = headerline(1:equal_pos-1)                          ;
            edfheader(noheaderlines,1) = {strtrim(left_str)}               ; %#ok<*AGROW>
            edfheader(noheaderlines,2) = {strtrim(right_str)}              ;
        end
    end
end
out.edfheader = edfheader                                                  ;
clear noheaderlines left_str right_str equal_pos semicolon_pos             ;
clear check_startofheader check_endofheader headerline edfheader           ;


% =========================================================================
% --- close fid for the moment to open later in the correct machine format
% =========================================================================
fclose(fid)                                                                ;


% =========================================================================
% --- determine some standard parameters from the edfheader
% =========================================================================
for i = 1 : size(out.edfheader,1)
    if ( strcmp(out.edfheader(i,1),'ByteOrder') == 1 )
        if ( strcmp(out.edfheader(i,2),'LowByteFirst') == 1 )
            HighByteFirst = 0                                              ;
        elseif ( strcmp(out.edfheader(i,2),'HighByteFirst') == 1 )
            HighByteFirst = 1                                              ;
        end
    end
    %---
    if ( strcmp(out.edfheader(i,1),'DataType') == 1 )
        if ( strcmp(out.edfheader(i,2),'UnsignedShort') == 1 )
            DataType = 0                                                   ;
            bytes    = 2                                                   ;
        elseif ( strcmp(out.edfheader(i,2),'UnsignedLong') == 1 )
            DataType = 1                                                   ;
            bytes    = 4                                                   ;
        elseif ( strcmp(out.edfheader(i,2),'SignedInteger') == 1 )
            DataType = 2                                                   ;
            bytes    = 4                                                   ;
        elseif ( strcmp(out.edfheader(i,2),'UnsignedByte') == 1 )          ;
            DataType = 4                                                   ;
            bytes    = 1                                                   ;
        end
    end
    %---
    if ( strcmp(out.edfheader(i,1),'Size') == 1 )
        dlen = str2num(out.edfheader{i,2})                                 ; %#ok<*ST2NM>
    end
    %---
    if ( strcmp(out.edfheader(i,1),'Dim_1') == 1 )
        cols = str2num(out.edfheader{i,2})                                 ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'Dim_2') == 1 )
        rows = str2num(out.edfheader{i,2})                                 ;
    end
end

% =========================================================================
% --- determine some extra parameters from the edfheader
% =========================================================================
for i = 1 : size(out.edfheader,1)
    if ( strcmp(out.edfheader(i,1),'col_beg') == 1 )
        col_beg = str2num(out.edfheader{i,2}) + 1                          ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'col_end') == 1 )
        col_end = str2num(out.edfheader{i,2}) + 1                          ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'row_beg') == 1 )
        row_beg = str2num(out.edfheader{i,2}) + 1                          ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'row_end') == 1 )
        row_end = str2num(out.edfheader{i,2}) + 1                          ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'count_time') == 1 )
        preset = str2num(out.edfheader{i,2})                               ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'prefix') == 1 )
        prefix = out.edfheader{i,2}                                        ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'suffix') == 1 )
        suffix = out.edfheader{i,2}                                        ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'run') == 1 )
        number = str2num(out.edfheader{i,2})                               ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'counter_mne') == 1 )
        counter_mne = out.edfheader{i,2}                                   ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'counter_pos') == 1 )
        counter_position = str2num(out.edfheader{i,2})                     ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'count_time') == 1 )
        count_time = str2num(out.edfheader{i,2})                           ;
    end
    %---
    if ( strcmp(out.edfheader(i,1),'time_of_frame') == 1 )
        time_of_frame = str2num(out.edfheader{i,2})                        ;
    end
end


% =========================================================================
% --- create variables, if certain expected header information were missing
% =========================================================================
if (  exist('number','var') == 0 )
    number = []                                                            ;
end
if (  exist('row_beg','var') == 0 )
    row_beg = 1                                                            ;
end
if (  exist('row_end','var') == 0 )
    row_end = rows                                                         ;
end
if (  exist('col_beg','var') == 0 )
    col_beg = 1                                                            ;
end
if (  exist('col_end','var') == 0 )
    col_end = cols                                                         ;
end
if (  exist('prefix','var') == 0 )
    prefix = ''                                                            ;
end
if (  exist('suffix','var') == 0 )
    suffix = '.edf'                                                        ;
end
if (  exist('preset','var') == 0 )
    preset = []                                                            ;
end


% =========================================================================
% --- try to find the elapsed time between frames
% =========================================================================
elapsed = []                                                               ;
if (  exist('counter_mne','var')==1 && exist('counter_position','var')==1  ...
   && ~isempty(counter_mne)   ==1 && ~isempty(counter_position)   ==1 )   
    counter_mne   = strtrim(counter_mne)                                   ;
    findSpace     = strfind(counter_mne,' ')                               ;
    counter_names = cell(numel(counter_position,1))                        ;
    for i = 1 : numel(counter_position)
        if ( i == 1 )
            counter_names{i,1} = strtrim(counter_mne(1:findSpace(1)-1))                    ;
        elseif ( i ~= 1 && i == numel(counter_position) )
            counter_names{i,1} = strtrim(counter_mne(findSpace(end)+1:numel(counter_mne))) ;
        else
            counter_names{i,1} = strtrim(counter_mne(findSpace(i-1)+1:findSpace(i)-1))     ;
        end
    end
    % ---
    for i = 1 : numel(counter_position)
        if ( strcmp(counter_names{i,1},'ccdt') == 1 )
            elapsed = counter_position(i)                                  ;
        end
    end
elseif  ( exist('count_time','var') == 1 && ~isempty(count_time) ==1 )
    elapsed =  count_time                                                  ; 
elseif  ( exist('time_of_frame','var') == 1 && ~isempty(time_of_frame) ==1)
    elapsed =  time_of_frame                                               ;
end


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
% --- Reopen the file in the right machine format and skip the header
% =========================================================================
if ( HighByteFirst == 0 )
    [fid,message] = fopen(file,'r')                                        ;            
elseif ( HighByteFirst == 1 )
    [fid,message] = fopen(file,'r','b')                                    ;            
end
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end

% --- Find the end of the header position
j = 1                                                                      ;
while j < 8
    seekpos  = 512*j - 2                                                   ; % the numbering starts at zero and a newline character follows
    status   = fseek(fid,seekpos,'bof')                                    ; %#ok<NASGU> % go to the next possible position
    newchar  = fread(fid,1,'char=>char')                                   ; % read the next Byte of the header   
    closecurlybracket = strfind(newchar,'}')                               ; % check if newchar matches '}'
    if ~isempty(closecurlybracket)                                           % if the close curly bracket '}' is reached break this loop
        break
    end
    j = j + 1                                                              ;
end
if ( j == 8 )
    uiwait(msgbox(message,'Could not find the end of the header!','error','modal')) ;
    return                                                                 ;
end
endofheaderpos = ftell(fid) + 1                                            ;

% --- Calculate the starting position of the data block from the file size
status   = fseek(fid,0,'eof')                                              ; %#ok<NASGU> % go to the end of the file
endpos   = ftell(fid)                                                      ; % get the block number of the end of the file
startpos = endpos - cols * rows * bytes                                    ; % estimate the start of the data section

% --- Compare the position values for the end of the header
if (startpos ~= endofheaderpos )
    uiwait(msgbox(message,'File seems to be corrupted!','error','modal'))  ;
    return                                                                 ;
end
clear j seekpos status newchar closecurlybracket endofheaderpos endpos     ;

% =========================================================================
% --- store image to output structure
% =========================================================================
status = fseek(fid,startpos,'bof')                                         ; %#ok<NASGU> % go to the start of the data block
if ( DataType == 0 )
    out.imm = fread(fid,[cols rows],'uint16=>single')                      ; % convert to single precison
elseif ( DataType == 1 )
    out.imm = fread(fid,[cols rows],'uint32=>single')                      ; % convert to single precison
elseif ( DataType == 2 )
    out.imm = fread(fid,[cols rows],'int32=>single')                       ; % convert to single precison
elseif ( DataType == 4 )
    out.imm = fread(fid,[cols rows],'uint8=>single')                       ; % convert to single precison
end
out.imm = out.imm'                                                         ; % the x and y coordinets are transposed
clear startpos                                                             ;

% =========================================================================
% --- close fid and return
% =========================================================================
fclose(fid)                                                                ;


% =========================================================================
% --- remove some variables
% =========================================================================
clear rows row_beg row_end cols col_beg col_end preset dlen suffix prefix  ;
clear DataType HighByteFirst bytes number i file fid message               ;
clear elapsed count_time varargin                                          ;


% ---
% EOF
