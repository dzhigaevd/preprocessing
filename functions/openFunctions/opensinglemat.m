function out = opensinglemat(varargin)
% ---
% --- OPENSINGLEMAT Open *.mat single files
% ---
% --- USAGE : 
% --- OUT = OPENSINGLEMAT(FILENAME)        --> direct call to a file
% --- OUT = OPENSINGLEMAT(FILENAME,INDEX)  --> call to open a specific 
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
% --- OUT : structure containing matheader, header & image. Usually image 
% ---       is an intensity image with 0 for black, saturation point for 
% ---       white (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2014/06/10 $ function to open *.mat files
% ---


% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
out = []                                                                   ;
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to a single MAT file
elseif ( nargin == 2 )
    firstbatchimage    = varargin{1}                                       ;
    Index              = varargin{2}                                       ;
    [pathstr,name,ext] = fileparts(firstbatchimage)                        ;
    findunderscore     = strfind(name,'_')                                 ;
    if ( ~isempty(findunderscore) == 1 )
        indexlastunderscore = findunderscore(end)                          ;
        namestring          = name(1:indexlastunderscore)                  ;
        namenumber          = sprintf('%05i',Index)                        ;
        file                = fullfile(pathstr,[namestring,namenumber,ext]); % create fullfile name of wanted image
    else
        uiwait(msgbox('MAT File Open Error','error','modal'))              ;
        return                                                             ;        
    end
end


% =========================================================================
% --- check if file exists
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
fclose(fid)                                                                ;

% =========================================================================
% --- load image & mat header from file
% --- structure 'out' with fields .imm' & '.matheader'
% =========================================================================
load(file)                                                                 ;

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
header(5,:)  = {'number',        1}                                        ;
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
out.header = header                                                        ;


% =========================================================================
% --- remove some variables
% =========================================================================
clear header rows cols bytes fid message varargin                          ;


% ---
% EOF
