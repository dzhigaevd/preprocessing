function out = opensinglecbf(varargin)
% ---
% --- OPENSINGLECBF Open *.cbf single files
% ---
% --- USAGE : 
% --- OUT = OPENSINGLECBF(FILE)        --> direct call to a file
% --- OUT = OPENSINGLECBF(FILE,INDEX)  --> call to open a specific 
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
% --- OUT : structure containing CBFheader, header & image. Usually image 
% ---       is an intensity image with 0 for black, saturation point for 
% ---       white (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2010/03/28 $ function to open *.cbf files
% ---
% --- Based on: 'cbfread.m' by the Swiss Light Source
% --- $Revision: 1.1 $  $Date: 2008/06/10 17:05:13 $
% --- $Author: bunk $
% --- $Tag: $
% ---
% --- Description:
% --- Macro for reading Crystallographic Binary File (CBF) files written
% --- by the Pilatus detector control program camserver. 
% ---
% --- Note:
% --- Compile the C-program cbf_uncompress using mex (see header of
% --- cbf_uncompress.c) to use it for uncompression instead of the slower
% --- Matlab code. 
% --- Currently this routine supports only the subset of CBF features
% --- needed to read the Pilatus detector data. 
% --- Call without arguments for a brief help text.
% ---
% --- Dependencies:
% ---  - image_read_set_default
% ---  - fopen_until_exists
% ---  - get_hdr_val
% ---  - compiling cbf_uncompress.c increases speed but is not mandatory
% ---

% =========================================================================
% --- Some additional flags and parameters
% =========================================================================
max_header_length = 4096                                                   ; % expected maximum length for the text header
eoh_signature     = char([ 12 26 4 213 ])                                  ; % end of header signature
cbf_signature     = '###CBF: VERSION'                                      ; % CBF file signature


% =========================================================================
% --- Check minimum number of input arguments
% =========================================================================
if (nargin < 1)
    uiwait(msgbox('CBF File Open Error','error','modal'))                  ;
    return                                                                 ;
end


% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to a single CBF file
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
        uiwait(msgbox('CBF File Open Error','error','modal'))              ;
        return                                                             ;
    end
end


% =========================================================================
% --- Calling an external C routine for uncompressing the data did save
% --- about 30% time on a specific machine. 
% --- The C-routine is used if a compiled version of it exists. 
% --- See the header of cbf_uncompress.c for information on how to compile
% --- the C file using mex in Matlab. 
% =========================================================================
c_routine = (length(which('cbf_uncompress')) > 0)                          ; %#ok<ISMT>


% =========================================================================
% --- open file, read all the data & close the file
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
fdat = fread(fid,'uint8=>uint8')                                           ;
fclose(fid)                                                                ;


% =========================================================================
% --- search for end of header signature within the expected maximum 
% --- length of a header
% =========================================================================
end_of_header_pos = strfind(fdat(1:min(max_header_length,length(fdat)))'...
                           ,eoh_signature )                                ;
if (length(end_of_header_pos) < 1)
    cbf_error(file,'no header end signature found')                        ;
    return                                                                 ;
end


% =========================================================================
% Call of subfunction: 'char_to_cellstr'
% Return the complete header as lines of a cell array
% Check for CBF signature
% =========================================================================
out.cbfheader = char_to_cellstr(char(fdat(1:(end_of_header_pos-1))'))      ;
if (~strncmp(cbf_signature,out.cbfheader{1},length(cbf_signature)))
    cbf_error(file,[ 'CBF signature ''' cbf_signature ''' not found in first line ''' out.cbfheader{1} '''' ]);
end


% =========================================================================
% Extract the mandatory information for decompression from the header
% =========================================================================
nbin_bytes = get_hdr_val(out.cbfheader,'X-Binary-Size:','%f',1)                   ;
dim1       = get_hdr_val(out.cbfheader,'X-Binary-Size-Fastest-Dimension:','%f',1) ;
dim2       = get_hdr_val(out.cbfheader,'X-Binary-Size-Second-Dimension:','%f',1)  ;
el_type    = get_hdr_val(out.cbfheader,'X-Binary-Element-Type: "','%[^"]',1)      ;
compr_type = get_hdr_val(out.cbfheader,'conversions="','%[^"]',1)                 ;
switch (el_type)
    case 'signed 32-bit integer'
        bytes = 4                                                          ;
    otherwise
        cbf_error(file,[ 'Unknown element type ' el_type ])                ;
end
switch (compr_type)
    case 'x-CBF_BYTE_OFFSET' 
        compression_type = 1                                               ;
    otherwise
        cbf_error(file,[ 'Unknown compression type ' compr_type ])         ;
end


% =========================================================================
% Extract usefull header information
% =========================================================================
elapsed = get_hdr_val(out.cbfheader,'# Exposure_period','%f',1)            ;
preset  = get_hdr_val(out.cbfheader,'# Exposure_time','%f',1)              ;
iname   = out.cbfheader{2,1}                                               ; % image name
findunderscore = strfind(iname,'_')                                        ;
if ( ~isempty(findunderscore) == 1 )
    indexlastunderscore = findunderscore(end)                              ;
    prefix    = iname(1:indexlastunderscore)                               ;
    number    = str2num(iname(indexlastunderscore+1:end))                  ; %#ok<ST2NM>
else
    uiwait(msgbox('CBF File Open Error: None standard name','error','modal'));
    return                                                                 ;
end
suffix = '.cbf'                                                            ;
clear findunderscore indexlastunderscore iname                             ;


% =========================================================================
% --- Uncompress the binary data
% =========================================================================
out.imm = extract_frame(fdat((end_of_header_pos+length(eoh_signature)):end)   ...
                       ,dim1,dim2,nbin_bytes,compression_type,file,c_routine)    ;
out.imm = flipud(out.imm')                                                 ; % correct image to natural view


% =========================================================================
% --- create field names of standard imm header structure
% --- leave values empty
% --- store header to output structure
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
header(9,:)  = {'row_beg',       1}                                        ;
header(10,:) = {'row_end',       size(out.imm,1)}                          ;
header(11,:) = {'col_beg',       1}                                        ;
header(12,:) = {'col_end',       size(out.imm,2)}                          ;
header(13,:) = {'row_bin',       1}                                        ;
header(14,:) = {'col_bin',       1}                                        ;
header(15,:) = {'rows',          size(out.imm,1)}                          ;
header(16,:) = {'cols',          size(out.imm,2)}                          ;
header(17,:) = {'bytes',         bytes}                                    ;
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
% ---
out.header = header                                                        ;
clear header                                                               ;


% =========================================================================
% --- Subfunction 'extract_frame'
% =========================================================================
function frame = extract_frame(dat_in,dim1,dim2,no_of_in_bytes          ...
                              ,compression_type,file,c_routine)

% =========================================================================
% --- If possible call external uncompression routine
% =========================================================================
if (c_routine)    
    frame=cbf_uncompress(dat_in,dim1,dim2,no_of_in_bytes,compression_type) ;
    return
end

% =========================================================================
% only byte-offset compression is supported
% =========================================================================
if (compression_type ~= 1)
    cbf_error(file,['extract_frame does not support compression type no. ' ...
                    num2str(compression_type)])                            ;
end

% =========================================================================
% --- In byte-offset compression the difference to the previous pixel value
% --- is stored as a byte, 16-bit integer or 32-bit integer, depending on 
% --- its size. 
% --- The sizes above one byte are indicated by the escape sequence -1 in 
% --- the previous data format, i.e, a 32-bit integer is preceded by the
% --- sequence 0x80 (too large for a byte) or 0x8000 (too large for a 
% --- 16-bit integer). 
% =========================================================================
frame    = zeros(dim1,dim2)                                                ;
ind_out  = 1                                                               ;
ind_in   = 1                                                               ;
val_curr = 0                                                               ;
val_diff = 0                                                               ; %#ok<NASGU>
while (ind_in <= no_of_in_bytes)
    val_diff = double(dat_in(ind_in))                                      ;
    ind_in = ind_in + 1                                                    ;
    if (val_diff ~= 128)
        % if not escaped as -128 (0x80=128) use the current byte as
        % difference, with manual complement to emulate the sign
        if (val_diff >= 128)
            val_diff = val_diff - 256                                      ;
        end
    else
        % otherwise check for 16-bit integer value
        if ((dat_in(ind_in) ~= 0) || (dat_in(ind_in+1) ~= 128))
            % if not escaped as -32768 (0x8000) use the current 16-bit integer as difference 
            val_diff = double(dat_in(ind_in))+256*double(dat_in(ind_in+1)) ;
            % manual complement to emulate the sign
            if (val_diff >= 32768)
                val_diff = val_diff - 65536                                ;
            end
            ind_in = ind_in + 2                                            ;
        else
            ind_in = ind_in + 2                                            ;
            % if everything else failed use the current 32-bit value as difference
            val_diff =            double(dat_in(ind_in))                   ...
                     +      256 * double(dat_in(ind_in+1))                 ...
                     +    65536 * double(dat_in(ind_in+2))                 ...
                     + 16777216 * double(dat_in(ind_in+3))                 ;
            % manual complement to emulate the sign
            if (val_diff >= 2147483648)
                val_diff = val_diff - 4294967296                           ;
            end
            ind_in = ind_in + 4                                            ;
        end
    end
	val_curr = val_curr + val_diff                                         ;
    frame(ind_out) = val_curr                                              ;
    ind_out = ind_out + 1                                                  ;
end
    
% =========================================================================
% --- Consistency check
% =========================================================================
if (ind_out-1 ~= dim1*dim2)
    cbf_error(file,['Mismatch between ' num2str(ind_out-1)                 ...
                    ' bytes after decompression with ' num2str(dim1*dim2)  ...
                    ' expected ones' ])                                    ;
end


% =========================================================================
% --- Subfunction 'cbf_error'
% =========================================================================
function [] = cbf_error(file,text)

fprintf('cbfread of %s:\n %s\n',file,text)                                 ;
return                                                                     ;


% =========================================================================
% --- Subfunction 'char_to_cellstr'
% =========================================================================
function [outstr] = char_to_cellstr(inchars,nl_only)

if (nargin < 2)
    nl_only = 0;
end

% =========================================================================
% --- get positions of end-of-line signatures
% =========================================================================
eol_ind  = regexp(inchars,'\r\n')                                          ;
eol_offs = 1                                                               ;
if ((length(eol_ind) < 1) || (nl_only))
    eol_ind  = regexp(inchars,'\n')                                        ;
    eol_offs = 0                                                           ;
end
if (length(eol_ind) < 1)
    eol_ind = length(inchars) + 1                                          ;
end
if (length(eol_ind) < 1)
    outstr = []                                                            ;
    return                                                                 ;
end

% =========================================================================
% dimension return array with number of lines
% =========================================================================
outstr = cell(length(eol_ind),1)                                           ;

% =========================================================================
% copy the lines to the return array, suppressing empty lines
% =========================================================================
start_pos = 1                                                              ;
ind_out   = 1                                                              ;
for ind = 1 : length(eol_ind)
    end_pos = eol_ind(ind) - 1                                             ;
    while ((end_pos >= start_pos) && (inchars(end_pos) == ' '))
        end_pos = end_pos - 1                                              ;
    end
    % store non-empty strings
    if (end_pos >= start_pos)
        outstr{ind_out} = inchars(start_pos:end_pos)                       ;
        ind_out = ind_out + 1                                              ;
    end
    start_pos = eol_ind(ind) + 1 + eol_offs                                ;
end

% =========================================================================
% resize cell array in case of empty lines
% =========================================================================
if (ind_out <= length(eol_ind))
    outstr = outstr(1:(ind_out-1))                                         ;
end

% =========================================================================
% --- Subfunction 'get_hdr_val'
% =========================================================================
% --- Description:
% --- Find text signature in a bunch of cell strings from a file header and
% --- return the following value in the specified format. 
% --- Example: 
% --- no_of_bin_bytes = get_hdr_val(header,'X-Binary-Size:','%f',1);
% --- The last parameter specifies whether the macro should exit with an 
% --- error message if the text signature has not been found. 
function [outval,line_number,err] = get_hdr_val(header,signature,format,exit_if_not_found)

% =========================================================================
% ---Initialize output arguments
% =========================================================================
outval      = 0                                                            ;
line_number = 0                                                            ;
err         = 0                                                            ;

% =========================================================================
% --- Search for the signature string
% =========================================================================
pos_found        = strfind(header,signature)                               ;
signature_sscanf = strrep(signature,'%','%%')                              ; % for sscanf the percentage sign has a special meaning

for ind=1:length(pos_found)                                                 % loop over the search results for all header lines
    if (length(pos_found{ind}) > 0)                                         %#ok<ISMT>
        % get the following value in the specified format
        [outval,count] = sscanf(header{ind}(pos_found{ind}:end),[signature_sscanf format]);
        if (count < 1)
            outval = 0                                                     ;
            err    = 1                                                     ;
        else
            if (count > 1)
                outval = outval(1)                                         ;
            end
            line_number = ind                                              ;
            return                                                         ;
        end
    end
end

% no occurrence found
err = 1                                                                    ;
if (exit_if_not_found)
    error(['No header line with signature ''' signature ''' and format ' format ' found']);
end
return


% ---
% EOF
