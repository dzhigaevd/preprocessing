function out = openh5attribute(varargin)

out = [];       

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
file               = varargin{1}                                           ;
attribute          = varargin{2}                                           ;

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

if master == 1
    dataInfo = h5info(file, attribute);
    
    dataFullSize  = dataInfo.Dataspace.Size;
    out = h5read(file,sprintf(attribute));
    
    %Specific for Ni6602_buff
    out = transpose(out);
    out = out(:,any(out));
    
    %To mimic the matrix structure in HDFView
%     switch dataFullSize
%         case 2
%             out = transpose(out); 
%         case {3,4}
%             out = permute(out,[3 1 2]); 
%     end
%     if attribute == 'Ni6602_buff'
%         out=out(:,any(out)) % Removes columns/rows with only zeros
%     end
end

