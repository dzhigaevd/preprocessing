function out = openmultispe(varargin)
% ---
% --- OPENMULTISPE Open *.spe multi files
% ---
% --- USAGE : 
% --- OUT = OPENMULTISPE(FILENAME,INDEX)  --> call to open a specific 
% --- image from a multi *.spe file. Here filename stands for the 
% --- full file name of the multifile and the index points to one of the  
% --- images of the multi *.spe file.
% ---
% --- Input Argument:
% --- FILENAME : image file name
% --- Index    : image index
% ---
% --- Output Argument:
% --- OUT : structure containing speheader, header & image. Usually image 
% ---       is an intensity image with 0 for black, saturation point for 
% ---       white (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2010/05/31 $ function to open multi *.spe files
% ---


% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to the first image of an multi SPE file
    Index    = 1                                                           ; % assume the first image is wanted
elseif ( nargin == 2 )
    file     = varargin{1}                                                 ;
    Index    = varargin{2}                                                 ;
end


% =========================================================================
% --- if image index is out of range, error and return
% =========================================================================
nimages = openspeheader(file)                                              ; % check how many images are contained in the multi file
if ( Index < 1 ||  Index > nimages )
    uiwait(msgbox('Image index is out of range!','OPENMULTISPE Error'   ...
                 ,'error','modal'))                                        ;
    return                                                                 ;
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
% --- load information from the first spe header
% =========================================================================
fseek(fid, 108,'bof'); DataType = fread(fid,1,'uint16=>single')            ;
fseek(fid,  42,'bof'); cols     = double(fread(fid,1,'uint16=>single'))    ; 
fseek(fid, 656,'bof'); rows     = double(fread(fid,1,'uint16=>single'))    ;


% =========================================================================
% --- determine header size & image size
% =========================================================================
headersize = 4100                                                          ; % spe standard header size
switch DataType
    case 0
        bytes   = 4                                                        ;
    case 1
        bytes   = 4                                                        ;
    case 2
        bytes   = 2                                                        ;
    case 3
        bytes   = 2                                                        ;
end
imagesize = bytes * rows *cols                                             ; % image size in bytes
clear bytes rows cols                                                      ;


% =========================================================================
% --- determine image start byte (only one header in multi *.spe files)
% =========================================================================
imagestartbyte  = (Index - 1) * imagesize + headersize                     ;


% =========================================================================
% --- load information from spe header
% =========================================================================
fseek(fid,   6,'bof'); speheader.xDimDet  = fread(fid,1,'uint16=>single')  ; 
fseek(fid,  10,'bof'); speheader.exposure = fread(fid,1,'float32=>single') ; % exposure time [s] 
fseek(fid,  18,'bof'); speheader.yDimDet  = fread(fid,1,'uint16=>single')  ; 
fseek(fid,  42,'bof'); speheader.xdim     = fread(fid,1,'uint16=>single')  ; 
fseek(fid, 108,'bof'); speheader.DataType = fread(fid,1,'uint16=>single')  ;
fseek(fid, 656,'bof'); speheader.ydim     = fread(fid,1,'uint16=>single')  ;
fseek(fid, 672,'bof'); speheader.readout  = fread(fid,1,'float32=>single') ; % readout time [ms]!!!  
fseek(fid,1446,'bof'); speheader.nimages  = fread(fid,1,'uint16=>single')  ; 
out.speheader = speheader                                                  ; % assign ouput variable


% =========================================================================
% --- load image from *.spe multi file
% =========================================================================
fseek(fid,imagestartbyte,'bof')                                            ; % jump to data start values
% ---	
DataType = speheader.DataType                                              ;
cols     = double(speheader.xdim)                                          ;        
rows     = double(speheader.ydim)                                          ;        
switch DataType
    case 0
        out.imm = fread(fid,[cols rows],'float32=>single')                 ; % FLOATING POINT (4 bytes / 32 bits)
        bytes   = 4                                                        ;
    case 1
        out.imm = fread(fid,[cols rows],'int32=>single')                   ; % LONG INTEGER (4 bytes / 32 bits)
        bytes   = 4                                                        ;
    case 2
        out.imm = fread(fid,[cols rows],'int16=>single')                   ; % INTEGER (2 bytes / 16 bits)
        bytes   = 2                                                        ;
    case 3
        out.imm = fread(fid,[cols rows],'uint16=>single')                  ; % UNSIGNED INTEGER (2 bytes / 16 bits)
        bytes   = 2                                                        ;
end
out.imm    = out.imm'                                                      ; % the x and y coordinets are transposed


% =========================================================================
% --- create field names of standard imm header structure
% --- leave values empty
% =========================================================================
header = cell(53,2)                                                        ; % initialize '.header' structure
header(1,:)  = {'mode',          []}                                       ;
header(2,:)  = {'compression',   []}                                       ; 
header(3,:)  = {'date',          []}                                       ;
header(4,:)  = {'prefix',        ''}                                       ;
header(5,:)  = {'number',        Index}                                    ;
header(6,:)  = {'suffix',        '.spe'}                                   ;
header(7,:)  = {'monitor',       []}                                       ;
header(8,:)  = {'shutter',       []}                                       ;
header(9,:)  = {'row_beg',       1}                                        ;
header(10,:) = {'row_end',       rows}                                     ;
header(11,:) = {'col_beg',       1}                                        ;
header(12,:) = {'col_end',       cols}                                     ;
header(13,:) = {'row_bin',       []}                                       ;
header(14,:) = {'col_bin',       []}                                       ;
header(15,:) = {'rows',          rows}                                     ;
header(16,:) = {'cols',          cols}                                     ;
header(17,:) = {'bytes',         bytes}                                    ;
header(18,:) = {'kinetics',      []}                                       ;
header(19,:) = {'kinwinsize',    []}                                       ;
header(20,:) = {'elapsed',       speheader.exposure+speheader.readout/1000};
header(21,:) = {'preset',        speheader.exposure}                       ;
header(22,:) = {'topup',         []}                                       ;
header(23,:) = {'inject',        []}                                       ;
header(24,:) = {'dlen',          cols * rows / bytes}                      ; % assume it is uncompressed
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
out.header = header                                                        ; % assign emulated imm header


% =========================================================================
% --- close fid and return
% =========================================================================
fclose(fid)                                                                ;


% =========================================================================
% --- remove some variables
% =========================================================================
clear header speheader rows cols bytes DataType fid message varargin       ;


% ---
% EOF
