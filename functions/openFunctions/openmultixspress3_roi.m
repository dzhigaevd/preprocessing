function out = openmultixspress3_roi(varargin)
% --- USAGE
% --- OUT = OPENMULTIEIGER4M(FILENAME)        -- direct call to 
% --- open all images of a HDF5 file
% --- OUT = OPENMULTIEIGER4M(FILENAME,STARTINDEX,ENDINDEX)-- call to 
% --- open a block of images (lines). 
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
% --- ROIextended [roiYstart, roiYend, roiXstart, roiXend, frameStart, frameEnd]

% --- Output Argument
% --- OUT  structure containing header & image block.
% ---
% --- D. Dzhigaev
% --- $Revision 1.0 $Date 20190411 $ function to open merlin files
% --- Adaptation from P10, PETRA III scripts by M. Sprung

out = [];       

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 2 )
    file               = varargin{1}                                       ;
    scan_meta          = varargin{2}                                       ;
    StartIndex         = 1                                                 ;
    EndIndex           = inf                                               ;
    roiYstart          = 1                                                 ;
    roiXstart          = 1                                                 ;
    roiYrange          = inf                                               ;
    roiXrange          = inf                                               ;
    frameStart         = 1                                                 ;
    frameEnd           = inf                                               ;
    
elseif ( nargin == 3 )
    file               = varargin{1}                                       ;
    scan_meta          = varargin{2}                                       ;
    StartIndex         = 1                                                 ;
    EndIndex           = inf                                               ;
    roiYstart          = varargin{3}(1)                                    ;
    roiYend            = varargin{3}(2)                                    ;
    roiXstart          = varargin{3}(3)                                    ;
    roiXend            = varargin{3}(4)                                    ;
    roiYrange          = roiYend - roiYstart + 1                           ;
    roiXrange          = roiXend - roiXstart + 1                           ;
    clear roiXend roiYend                                                  ;
elseif ( nargin == 4 )
    file               = varargin{1}                                       ;
    scan_meta          = varargin{2}                                       ;
    StartIndex         = varargin{3}                                       ;
    EndIndex           = varargin{4}                                       ;
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
    
elseif ( nargin == 5 )
    file               = varargin{1}                                       ;
    scan_meta          = varargin{2}                                       ;
    StartIndex         = varargin{3}                                       ;
    EndIndex           = varargin{4}                                       ;
    if StartIndex >= EndIndex
        dummy      = EndIndex                                              ;
        EndIndex   = StartIndex                                            ;
        StartIndex = dummy                                                 ;
        clear dummy                                                        ;
    end
    roiYstart          = varargin{5}(1)                                    ;
    roiYend            = varargin{5}(2)                                    ;
    roiXstart          = varargin{5}(3)                                    ;
    roiXend            = varargin{5}(4)                                    ;
    frameStart         = varargin{5}(5)                                    ;
    frameEnd           = varargin{5}(6)                                    ;
    
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
if strcmp(ext,'.h5')
    master     = 1                                                         ;
    masterfile = file                                                      ;
    disp('Using master file!');
elseif strcmp(ext,'.hdf5')
    data = 1                                                               ;
    datafile = file;
    disp('Using data file!');
else
    disp('Error! Input File not recognized as master or data file')        ;
    return
end

if data == 1
    fileInfo = h5info(file)                                                ;
    dataInfo = h5info(file, '/entry_0000/measurement/xspress3/data/')        ;
    
    % dataFullSize: [xDiff,yDiff,fastAxis,slowAxis]
    dataFullSize  = [dataInfo.Dataspace.Size, length(fileInfo.Groups)]     ;
    
    if EndIndex == inf
        EndIndex = dataFullSize(4)                              ; % set EndIndex to the chunk size
    end
    
    if frameEnd == inf
        frameEnd = dataFullSize(3)                              ; % set EndIndex to the chunk size
    end
    
    if dataInfo.Dataspace.Size(3) >= StartIndex %&& dataInfo.Dataspace.Size(3) >= EndIndex                    
        try
            switch scan_meta.type
                case 'flyscan'
                    out = zeros([dataFullSize(1),dataFullSize(2),frameEnd - frameStart+1,EndIndex - StartIndex+1],'single');
                case 'dmesh'
                    out = zeros([dataFullSize(1),dataFullSize(2),scan_meta.dim0,scan_meta.dim1],'single');
            end
        catch EM
            disp(EM)                                               ;
            fprintf('Please load the data in smaller chunks! \n')     ;
            return                                                 ;
        end
        
        for ii = 1 : EndIndex - StartIndex + 1
            for jj = 1 : frameEnd - frameStart + 1                
                out(:,:,jj,ii) = h5read(file,sprintf('/entry_%04i/measurement/xspress3/data/',StartIndex+ii-2),[roiXstart roiYstart frameStart+jj-1],[roiXrange roiYrange 1]);                
            end
            fprintf('Read line #%d \n',ii);
        end
        
        if strcmp(scan_meta.type,'dmesh')
            out_old = out;
            clear out;
            counter = 1;
            for i=1:scan_meta.dim1 %21
                for j=1:scan_meta.dim0 %41
                    out(:,:,j,i) = out_old(:,:,1,counter);
                    %fprintf('Moving datapoint %d to [%d,%d] \n',counter,j,i);
                    counter = counter + 1;
                end
            end
        end
        
        fprintf('Final data size [%d,%d,%d,%d] \n', size(out));
    else
        fprintf('Error! Data file contains only %i frames!',dataInfo.Dataspace.Size(3)) ; %DSPS
        return
    end 
end


end

