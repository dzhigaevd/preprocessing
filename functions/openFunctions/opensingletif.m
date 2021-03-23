function out = opensingletif(varargin)
% ---
% --- OPENSINGLETIF Open *.tif single files
% ---
% --- USAGE : 
% --- OUT = OPENSINGLETIF(FILENAME)        --> direct call to a file
% --- OUT = OPENSINGLETIF(FILENAME,INDEX)  --> call to open a specific 
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
% --- OUT : structure containing tifheader, header & image. Usually image 
% ---       is an intensity image with 0 for black, saturation point for 
% ---       white (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2008/08/05 $ function to open either compressed
% ---                                    or uncompresses *.tif single files
% ---


% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to a single TIF file
elseif ( nargin == 2 )
    firstbatchimage    = varargin{1}                                       ;
    Index              = varargin{2}                                       ;
    [pathstr,name,ext] = fileparts(firstbatchimage)                        ;
    namestring         = name(1:end-5)                                     ; % remove the number part behind the last underscore
    namenumber         = sprintf('%05i',Index)                             ; % create new 4 digit number part of the name string
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
% --- read tif data header (4096 Bytes)
% =========================================================================
tifheader = fread(fid,4096,'*char')'                                       ;
% ---
pound_pos = find(tifheader == '#')                                         ;
if ( ~isempty(pound_pos) == 1 )
    header_info = cell(numel(pound_pos),1)                                 ;
    for i = 1 : numel(pound_pos)-1
        header_info{i} = tifheader(pound_pos(i):pound_pos(i+1)-1)          ;
    end
    last_pound = tifheader(pound_pos(end):end)                             ;
    header_info{numel(pound_pos)} = char(cellstr(last_pound))              ;
    % ---
    for i = 1 : numel(pound_pos)
        new_str        = strtrim(header_info{i})                           ;
        header_info{i} = new_str(3:end)                                    ;
    end
    % ---
    clear last_pound pound_pos new_str                                     ;
    % --- get some header information
    for i = 1 : numel(header_info)
        if ( strcmp('Exposure_time',header_info{i}) == 1 )
            space_pos = find(header_info{i} == ' ')                        ;
            preset    = str2num(header_info{i}(space_pos(1):space_pos(2))) ;  %#ok<ST2NM>
        elseif ( strcmp('Exposure_period',header_info{i}) == 1 )
            space_pos = find(header_info{i} == ' ')                        ;
            elapsed   = str2num(header_info{i}(space_pos(1):space_pos(2))) ;  %#ok<ST2NM>
        end
    end
    % ---
end
% ---
out.tifheader = header_info                                                ;
clear header_info                                                          ;


% =========================================================================
% --- import tif data
% =========================================================================
out.imm = single(importdata(file))                                         ;
out.imm = out.imm'                                                         ; % use if the x and y coordinets are transposed


% =========================================================================
% --- close fid
% =========================================================================
fclose(fid)                                                                ;


% =========================================================================
% --- create field names of standard imm header structure
% --- leave unknown values empty
% =========================================================================
header = cell(53,2)                                                        ; % initialize '.header' structure
header(1,:)  = {'mode',          []}                                       ;
header(2,:)  = {'compression',   []}                                       ; 
header(3,:)  = {'date',          []}                                       ;
header(4,:)  = {'prefix',        []}                                       ;
header(5,:)  = {'number',        []}                                       ;
header(6,:)  = {'suffix',        []}                                       ;
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
header(24,:) = {'dlen',          []}                                       ; % assume it is uncompressed
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
clear header preset elapsed                                                ;


% ---
% EOF
