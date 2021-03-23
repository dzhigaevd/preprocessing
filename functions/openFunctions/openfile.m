function out = openfile(varargin)
% ---
% --- OPENFILE File opening function for different data file formats
% ---
% --- This function uses the file extension to separate between different
% --- file formats.
% ---
% --- OUT = OPENFILE(FILENAME,INDEX,{IMAGESTARTBYTE})
% ---
% --- DATA FILE NAMING CONVENTION :
% --- File names of single images CANNOT contain 'dashes'. The file name
% --- index should be written after an 'underscore' directly before the
% --- 'dot'.
% --- File names of multi image files should use a 'dash' to separate the
% --- first index number from the last index number. The first index
% --- number should be written after the last 'underscore' directly before
% --- the 'dot'.
% ---
% --- Input Argument(s):
% --- 1) FILENAME      : image file name
% ---
% --- Mandatory for multifiles:
% --- 2) INDEX         : image index number
% ---
% --- Optinal          :
% --- a) for compressed '*.imm' multifiles
% --- IMAGESTARTBYTE   : previously by INDEXCOMPRESSEDMULTIIMM determined
% ---                    starting byte positition of image #index
% --- b) to be continued
% ---
% ---
% --- Output Argument:
% --- out : structure containing header and image
% ---
% --- out.header : cell structure containing the header name & header value
% --- out.imm    : contains the image in single format
% ---
% ---  Michael Sprung
% ---  $Revision: 1.0 $  $Date: 2006/09/06 $
% ---


% =========================================================================
% --- prepare file name
% =========================================================================
file               = varargin{1}                                           ;
[pathstr,name,ext] = fileparts(file)                                       ; %#ok<ASGLU>


% =========================================================================
% --- check file mode (single or multi file)
% --- not perfect --> better to read first header and find switch
% =========================================================================
multi = 0                                                                  ; % assume single file format
% ---
switch (ext)
    case {'.imm', '.IMM'}
        findDash = strfind(name,'-')                                       ;
        if ( ~isempty(findDash) == 1 )                                       % file name does contain dashes --> multi images
            multi = 1                                                      ;
        end
    % ---    
    case {'.spe', '.SPE'}
         nimages = openspeheader(file)                                     ;
         if ( nimages > 1 )
            multi = 1                                                      ;             
         end
    % ---    
    case {'.h5', '.H5'}
          multi = 1                                                        ; % always multi for hdf5 files           
    % ---    
    case {'.nxs', '.NXS'}
          multi = 1                                                        ; % always multi for NEXUS (nxs) files           
    % ---
    case {'.mat', '.MAT'}
        findDash = strfind(name,'-')                                       ;
        if ( ~isempty(findDash) == 1 )                                       % file name does contain dashes --> multi images
            multi = 1                                                      ;
        end
    % ---
end


% =========================================================================
% --- check if input provides the (mandatory) Index for multifiles 
% =========================================================================
Index = []                                                                 ;
if ( multi == 1 )
    switch (ext)
        % ---
        case {'.imm', '.IMM'}
            if ( nargin > 1 )
                Index = varargin{2}                                        ;
            else
                Index = str2num(cell2mat(inputdlg(                      ...
                        'Please provide index number of image to open'  ...
                       ,'Openfile dialog')))                               ; 
            end
        % ---
        case {'.spe', '.SPE'}
            if ( nargin > 1 )
                Index = varargin{2}                                        ;
                if ( Index > nimages )
                    prompt = sprintf('Error: Multi *.spe file contains only %i images. Please provide a new index number!',nimages);
                    Index  = str2num(cell2mat(inputdlg(prompt,'Openfile dialog'))) ; 
                    clear prompt                                           ;
                end
            else
                Index = str2num(cell2mat(inputdlg(                      ...
                        'Please provide index number of image to open'  ...
                       ,'Openfile dialog')))                               ; 
            end
        % ---
        case {'.mat', '.MAT'}
            if ( nargin > 1 )
                Index = varargin{2}                                        ;
            else
                Index = str2num(cell2mat(inputdlg(                      ...
                        'Please provide index number of image to open'  ...
                       ,'Openfile dialog')))                               ; 
            end
        % ---
        case {'.h5', '.H5'}
            if ( nargin > 1 )
                Index = varargin{2}                                        ;
            else
                Index = str2num(cell2mat(inputdlg(                      ...
                        'Please provide index number of image to open'  ...
                       ,'Openfile dialog')))                               ; %#ok<*ST2NM>
            end
        % ---
        case {'.nxs', '.NXS'}
            if ( nargin > 1 )
                Index = varargin{2}                                        ;
            else
                Index = str2num(cell2mat(inputdlg(                      ...
                        'Please provide index number of image to open'  ...
                       ,'Openfile dialog')))                               ;
            end
        % ---
    end
else
    if ( nargin > 1 )
        Index = varargin{2}                                                ;
    end
end


% =========================================================================
% --- call the right open routine for the input file
% =========================================================================
if ( multi == 0 )                                                           % single file case --> only file name needed
    % ---
    switch (ext)
        case {'.imm', '.IMM'}
            if ( nargin == 1 )
                out = opensingleimm(file)                                  ; % for a direct call to a single IMM file
            elseif ( nargin == 2 )
                out = opensingleimm(file,Index)                            ; % for a call to the index single IMM file of a batch
            end
        case {'.mat', '.MAT'}
            if ( nargin == 1 )
                out = opensinglemat(file)                                  ; % for a direct call to a single mat file
            elseif ( nargin == 2 )
                out = opensinglemat(file,Index)                            ; % for a call to the index single mat file of a batch
            end
        case {'.edf', '.EDF'}
            if ( nargin == 1 )
                out = opensingleedf(file)                                  ; % for a direct call to a single EDF file
            elseif ( nargin == 2 )
                out = opensingleedf(file,Index)                            ; % for a call to the index single EDF file of a batch
            end
        case {'.spe', '.SPE'}
            if ( nargin == 1 )
                out = opensinglespe(file)                                  ; % for a direct call to a single SPE file
            elseif ( nargin == 2 )
                out = opensinglespe(file,Index)                            ; % for a call to the index single SPE file of a batch
            end
        case {'.cbf', '.CBF'}
            if ( nargin == 1 )
                out = opensinglecbf(file)                                  ; % for a direct call to a single CBF file
            elseif ( nargin == 2 )
                out = opensinglecbf(file,Index)                            ; % for a call to the index single CBF file of a batch
            end
        case {'.img', '.IMG'}
            if ( nargin == 1 )
                out = opensingleimg(file)                                  ; % for a direct call to a single IMG file
            elseif ( nargin == 2 )
                out = opensingleimg(file,Index)                            ; % for a call to the index single IMG file of a batch
            end
        case {'.tif','.tiff', '.TIF', '.TIFF'}
            if ( nargin == 1 )
                out = opensingletif(file)                                  ; % for a direct call to a single TIF file
            elseif ( nargin == 2 )
                out = opensingletif(file,Index)                            ; % for a call to the index single TIF file of a batch
            end
        case {'.txt','.TXT'}
            if ( nargin == 1 )
                out = opensingletxt(file)                                  ; % for a direct call to a single TXT file
            elseif ( nargin == 2 )
                out = opensingletxt(file,Index)                            ; % for a call to the index single TXT file of a batch
            end
    end
    % ---
elseif ( multi == 1 )                                                        % multifile case --> file name & image index needed
    % ---
    switch (ext)
        % ---
        case {'.imm', '.IMM'}
            if ( nargin <= 2 )
                out = openmultiimm(file,Index)                             ; % for a multi IMM file
            elseif ( nargin == 3 )
                ImageStartByte = varargin{3}                               ;
                out = openmultiimm(file,Index,ImageStartByte)              ; % for a multi IMM file                
            end
        % ---    
        case {'.spe', '.SPE'}
            if ( nargin <= 2 )
                out = openmultispe(file,Index)                             ; % for a multi SPE file
            end
        % ---    
        case {'.mat', '.MAT'}
            if ( nargin <= 2 )
                out = openmultimat(file,Index)                             ; % for a multi mat file
            elseif ( nargin == 3 )
                ImageStartByte = varargin{3}                               ;
                out = openmultimat(file,Index,ImageStartByte)              ; % for a multi mat file                
            end
        % ---
        case {'.h5', '.H5'}
            if ( nargin <= 2 )
                out    = openh5(file,Index)                                ; % for a h5 file
            end
        % ---
        case {'.nxs', '.NXS'}
            if ( nargin <= 2 )
                out    = opennxs(file,Index)                               ; % for a nxs file
            end
        % ---     
    end
    % ---
end


% ---
% EOF
