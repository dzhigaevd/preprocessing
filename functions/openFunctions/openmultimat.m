function out = openmultimat(varargin)
% ---
% --- OPENMULTIMAT Open *.mat multi files
% ---
% --- USAGE : 
% --- OUT = OPENMULTIMAT(FILENAME,INDEX)  --> call to open a specific 
% --- image from a multi *.mat file. Here filename stands for the 
% --- full file name of the multifile and the index points to one of the  
% --- images of the multi *.mat file.
% ---
% --- The convention is that multi *.mat files contain the variables
% --- tensor(rows, cols, nframes) and the cell tensorheader. This cell
% --- structure contains atleast the fields 'bytes', 'exposure' and
% --- 'elapsed'. The multi *.mat file is stored with the option '-v7.3'!
% --- e.g. with "save(filename,'tensor','tensorheader','-v7.3')"
% ---
% --- Input Argument:
% --- FILENAME : image file name
% --- Index    : image index
% ---
% --- Output Argument:
% --- OUT : structure containing matheader, header & imm. Usually image 
% ---       is an intensity image with 0 for black, saturation point for 
% ---       white (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2014/06/10 $ function to open multi *.mat files
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
% --- input file preparation
% =========================================================================
findUnderscore = strfind(file,'_')                                         ;
findDash       = strfind(file,'-')                                         ;
findDot        = strfind(file,'.')                                         ;
firstIndex     = str2num(file(findUnderscore(end)+1:findDash(end)-1))      ; %#ok<*ST2NM>
lastIndex      = str2num(file(findDash(end)+1:findDot(end)-1))             ;


% =========================================================================
% --- if image index is out of range, error and return
% =========================================================================
if ( Index < firstIndex || Index > lastIndex )
    error('Image index is out of range.')                                  ;
end


% =========================================================================
% --- correct image index basing on the total image numbers in mat file
% =========================================================================
Index = Index - firstIndex + 1                                             ;


% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
fclose(fid)                                                                ;

% =========================================================================
% --- load image & header from the multi mat file
% --- and create structure 'out' with fields '.imm' & '.matheader'
% =========================================================================
m             = matfile(file)                                              ;
out.imm(:,:)  = m.data(:,:,Index)                                          ;
header        = m.header                                                   ;
out.matheader = header{Index}                                              ;
clear m tensorheader

% =========================================================================
% --- assign variables from '*.matheader' & data image
% =========================================================================
rows    = size(out.imm,1)                                                  ;
cols    = size(out.imm,2)                                                  ;
bytes   = out.matheader.bytes                                              ; % used bytes per pixel
preset  = out.matheader.exposure                                           ; % exposure time [s]
elapsed = out.matheader.elapsed                                            ; % elapsed time since the exposure of the 1st image started

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
header(6,:)  = {'suffix',        '.mat'}                                   ;
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
header(20,:) = {'elapsed',       elapsed}                                  ;
header(21,:) = {'preset',        preset}                                   ;
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
% --- remove some variables
% =========================================================================
clear header rows cols bytes fid message varargin                          ;


% ---
% EOF
