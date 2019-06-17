function [data, counters, cstr] = readfio(filename)
% --- Function to read fio files
% ---
% --- Input : Full file name of the *.fio file
% --- Output: 
% ---            data     : matrix of data values
% ---            counters : cell of counter names
% ---

displayon   = 0                                                            ;
data        = []                                                           ;
counters    = cell(100,1)                                                  ;
cstr        = ''                                                           ;

fid         = fopen(filename,'r')                                          ; %#ok<*ST2NM>
counternum  = 0                                                            ;
linecounter = 1                                                            ;
while ~feof(fid)
    tline          = fgetl(fid)                                            ;
    [token,remain] = strtok(tline)                                         ;
    if (linecounter == 5)
        cstr = tline                                                       ;
        if displayon == 1
            fprintf('%s collected data for this command: %s\r'          ...
                   , filename, cstr)                                       ;
        end
    end
    if ( strcmp(token,'Col') == 1 )
        spaces = strfind(remain,' ')                                       ;
        if ( ~isempty(spaces) && numel(spaces) >=3 )
            counternum = counternum + 1                                    ;
            counters{counternum} = remain(spaces(2)+1:spaces(3)-1)         ;
        end
    end
    if ~isempty(str2num(token)) 
        currow = str2num(tline)                                            ;
        data   = [ data; currow]                                           ;
    end
    linecounter = linecounter + 1                                          ;
end
counters(counternum+1:end) = []                                            ;
fclose(fid)                                                                ;

% ---
% EOF
