function out = opensingleimm(varargin)
% ---
% --- OPENSINGLEIMM Open *.imm single files
% ---
% --- USAGE : 
% --- OUT = OPENSINGLEIMM(FILENAME)        --> direct call to a file
% --- OUT = OPENSINGLEIMM(FILENAME,INDEX)  --> call to open a specific 
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
% --- OUT : structure containing header & image. Usually image is an
% ---       intensity image with 0 for black, saturation point for white
% ---       (65535 for 16 bit images).
% ---
% --- Michael Sprung
% --- $Revision: 1.0 $Date: 2006/09/06 $ function to open either compressed
% ---                                    or uncompresses *.imm single files
% ---


% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file     = varargin{1}                                                 ; % for a direct call to a single IMM file
elseif ( nargin == 2 )
    firstbatchimage    = varargin{1}                                       ;
    Index              = varargin{2}                                       ;
    [pathstr,name,ext] = fileparts(firstbatchimage)                        ;
    findUnderscore     = strfind(name,'_')                                 ;
    namestring         = name(1:findUnderscore(end))                       ; % remove the number part behind the last underscore
    namenumber         = sprintf('%05i',Index)                             ; % create new 5 digit number part of the name string
    file               = fullfile(pathstr,[namestring,namenumber,ext])     ; % create fullfile name of wanted image
end


% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file)                                                ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end


% =========================================================================
% --- load header
% =========================================================================
fseek(fid,0,'bof')                                                         ;
header = cell(53,2)                                                        ; % initialize '.header' structure
header(1,:)  = {'mode',          fread(fid,1 ,'int')}                      ; % starting byte   0 : length  4 bytes
header(2,:)  = {'compression',   fread(fid,1 ,'int')}                      ; % starting byte   4 : length  4 bytes
header(3,:)  = {'date',          fread(fid,32,'*char')'}                   ; % starting byte   8 : length 32 bytes
header(4,:)  = {'prefix',        fread(fid,16,'*char')'}                   ; % starting byte  40 : length 16 bytes
header(5,:)  = {'number',        fread(fid,1 ,'int')}                      ; % starting byte  56 : length  4 bytes
header(6,:)  = {'suffix',        fread(fid,16,'*char')'}                   ; % starting byte  60 : length 16 bytes
header(7,:)  = {'monitor',       fread(fid,1 ,'int')}                      ; % starting byte  76 : length  4 bytes
header(8,:)  = {'shutter',       fread(fid,1 ,'int')}                      ; % starting byte  80 : length  4 bytes
header(9,:)  = {'row_beg',       fread(fid,1 ,'int')}                      ; % starting byte  84 : length  4 bytes
header(10,:) = {'row_end',       fread(fid,1 ,'int')-1}                    ; % starting byte  88 : length  4 bytes
header(11,:) = {'col_beg',       fread(fid,1 ,'int')}                      ; % starting byte  92 : length  4 bytes
header(12,:) = {'col_end',       fread(fid,1 ,'int')-1}                    ; % starting byte  96 : length  4 bytes
header(13,:) = {'row_bin',       fread(fid,1 ,'int')}                      ; % starting byte 100 : length  4 bytes
header(14,:) = {'col_bin',       fread(fid,1 ,'int')}                      ; % starting byte 104 : length  4 bytes
header(15,:) = {'rows',          fread(fid,1 ,'int')}                      ; % starting byte 108 : length  4 bytes
header(16,:) = {'cols',          fread(fid,1 ,'int')}                      ; % starting byte 112 : length  4 bytes
header(17,:) = {'bytes',         fread(fid,1 ,'int')}                      ; % starting byte 116 : length  4 bytes
header(18,:) = {'kinetics',      fread(fid,1 ,'int')}                      ; % starting byte 120 : length  4 bytes
header(19,:) = {'kinwinsize',    fread(fid,1 ,'int')}                      ; % starting byte 124 : length  4 bytes
header(20,:) = {'elapsed',       fread(fid,1 ,'double')}                   ; % starting byte 128 : length  8 bytes
header(21,:) = {'preset',        fread(fid,1 ,'double')}                   ; % starting byte 136 : length  8 bytes
header(22,:) = {'topup',         fread(fid,1 ,'int')}                      ; % starting byte 144 : length  4 bytes
header(23,:) = {'inject',        fread(fid,1 ,'int')}                      ; % starting byte 148 : length  4 bytes
header(24,:) = {'dlen',          fread(fid,1 ,'int')}                      ; % starting byte 152 : length  4 bytes
header(25,:) = {'roi_number',    fread(fid,1 ,'int')}                      ; % starting byte 156 : length  4 bytes
% --- modified position as of 20060306
header(26,:) = {'buffer_number', fread(fid,1 ,'uint32')}                   ; % starting byte 160 : length  4 bytes
header(27,:) = {'systick',       fread(fid,1 ,'uint32')}                   ; % starting byte 164 : length  4 bytes
% --- shifted header positions as of 20060306
header(28,:) = {'pv1',           fread(fid,40,'*char')}                    ; % starting byte 168 : length 40 bytes
header(29,:) = {'pv1VAL',        fread(fid,1 ,'float')}                    ; % starting byte 208 : length  4 bytes
header(30,:) = {'pv2',           fread(fid,40,'*char')}                    ; % starting byte 212 : length 40 bytes
header(31,:) = {'pv2VAL',        fread(fid,1 ,'float')}                    ; % starting byte 252 : length  4 bytes
header(32,:) = {'pv3',           fread(fid,40,'*char')}                    ; % starting byte 256 : length 40 bytes
header(33,:) = {'pv3VAL',        fread(fid,1 ,'float')}                    ; % starting byte 296 : length  4 bytes
header(34,:) = {'pv4',           fread(fid,40,'*char')}                    ; % starting byte 300 : length 40 bytes
header(35,:) = {'pv4VAL',        fread(fid,1 ,'float')}                    ; % starting byte 340 : length  4 bytes
header(36,:) = {'pv5',           fread(fid,40,'*char')}                    ; % starting byte 344 : length 40 bytes
header(37,:) = {'pv5VAL',        fread(fid,1 ,'float')}                    ; % starting byte 384 : length  4 bytes
header(38,:) = {'pv6',           fread(fid,40,'*char')}                    ; % starting byte 388 : length 40 bytes
header(39,:) = {'pv6VAL',        fread(fid,1 ,'float')}                    ; % starting byte 428 : length  4 bytes
header(40,:) = {'pv7',           fread(fid,40,'*char')}                    ; % starting byte 432 : length 40 bytes
header(41,:) = {'pv7VAL',        fread(fid,1 ,'float')}                    ; % starting byte 472 : length  4 bytes
header(42,:) = {'pv8',           fread(fid,40,'*char')}                    ; % starting byte 476 : length 40 bytes
header(43,:) = {'pv8VAL',        fread(fid,1 ,'float')}                    ; % starting byte 516 : length  4 bytes
header(44,:) = {'pv9',           fread(fid,40,'*char')}                    ; % starting byte 520 : length 40 bytes
header(45,:) = {'pv9VAL',        fread(fid,1 ,'float')}                    ; % starting byte 560 : length  4 bytes
header(46,:) = {'pv10',          fread(fid,40,'*char')}                    ; % starting byte 564 : length 40 bytes
header(47,:) = {'pv10VAL',       fread(fid,1 ,'float')}                    ; % starting byte 604 : length  4 bytes
% --- new fields (added 10/2006)
header(48,:) = {'imageserver',   fread(fid,1 ,'float')}                    ; % starting byte 608 : length  4 bytes
header(49,:) = {'CPUspeed',      fread(fid,1 ,'float')}                    ; % starting byte 612 : length  4 bytes
header(50,:) = {'immversion',    fread(fid,1 ,'int')}                      ; % starting byte 616 : length  4 bytes
header(51,:) = {'corecotick',    fread(fid,1 ,'uint32')}                   ; % starting byte 620 : length  4 bytes
header(52,:) = {'cameratype',    fread(fid,1 ,'int')}                      ; % starting byte 624 : length  4 bytes
header(53,:) = {'threshhold',    fread(fid,1 ,'float')}                    ; % starting byte 628 : length  4 bytes


% =========================================================================
% --- store header to output structure
% =========================================================================
out.header = header                                                        ;


% =========================================================================
% --- check if data is compressed
% =========================================================================
compression = 0                                                            ; % assume uncompressed files
if ( header{2,2} == 65540 && isempty(strfind(file,'_ucp_')) )                % condition for SMD legacy code
    compression = 1                                                        ;
elseif ( header{2,2} == 6 && header{50,2} >= 11 )                            % condition for compressed data from Imageserver program
    compression = 1                                                        ;
end


% =========================================================================
% --- store image to output structure
% =========================================================================
if ( compression == 0 )

    fseek(fid,1024,'bof')                                                  ;
    if ( header{50,2} >= 11 )
        if ( header{17,2} == 2 )
            out.imm = fread(fid,[header{16,2},header{15,2}],'uint16=>single')' ;
        elseif ( header{17,2} == 4 )
            out.imm = fread(fid,[header{16,2},header{15,2}],'uint32=>single')' ;
        end
    else
        out.imm = fread(fid,[header{16,2},header{15,2}],'uint16=>single')' ;
    end

elseif ( compression == 1 )

    if ( header{2,2} == 65540 && isempty(strfind(file,'_ucp_')) )            % condition for SMD legacy code
        PixelNumber = floor(double(header{24,2}/8))                        ; % for compressed mode 'dlen' is the data length in bytes
        if ( ispc == 1 )
            fseek(fid ,1024,'bof')                                         ; % go to the data start of the compressed file
            data       = fread(fid,2*PixelNumber,'uint32')                 ; % read the full data block 
            datavalues = data(2:2:end)                                     ; % create vector with intensity values
            PixelIndex = data(1:2:end-1) + 1                               ; % create vector with pixel indices
            PixelValue = single(uread(uwrite(datavalues,'uint32')       ...
                                     ,PixelNumber,'float32'))'             ; % convert data values to binary to float32 to single
        else
            fseek(fid ,1024,'bof')                                         ; % go to the data indices start of the compressed file
            PixelIndex = fread(fid,PixelNumber,'uint32',4) + 1             ; % read all pixel indices            
            fseek(fid ,1024 + 4,'bof')                                     ; % go to the data values start of the compressed file
            PixelValue = fread(fid,PixelNumber,'float32=>single',4)        ; % read all pixel values                        
        end
        % --- create uncompressed data
        UCData = zeros(header{16,2},header{15,2},'single')                 ;
        UCData(PixelIndex) = PixelValue                                    ;
        % --- store image data
        out.imm = UCData'                                                  ;
    elseif ( header{2,2} == 6 && header{50,2} >= 11 )
        PixelNumber = floor(double(header{24,2}))                          ; % for compressed mode 'dlen' is the # of exposed pixels
        fseek(fid ,1024,'bof')                                             ; % go to next data start of the compressed file
        PixelIndex = fread(fid,PixelNumber,'uint32') + 1                   ; % read the full index block (correction to start at 1)
        if ( header{17,2} == 2 )
            PixelValue = fread(fid,PixelNumber,'uint16=>single')           ; % read the full value block
        elseif ( header{17,2} == 4 )
            PixelValue = fread(fid,PixelNumber,'uint32=>single')           ; % read the full value block            
        end
        % --- create uncompressed data
        UCData = zeros(header{16,2},header{15,2},'single')                 ;
        UCData(PixelIndex) = PixelValue                                    ;
        % --- store image data
        out.imm = UCData'                                                  ;
    end
    
end


% =========================================================================
% --- close fid and return
% =========================================================================
fclose(fid)                                                                ;


% ---
% EOF
