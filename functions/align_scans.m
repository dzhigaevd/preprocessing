function cumulative_offset = align_scans(template, compare)
%ALIGN_SCANS Summary of this function goes here
%   Detailed explanation goes here

% cumulative_offset, a measurement of how good the "compare" fits to the
% template
% template is what scan we adjust after, accumulation of abs of offsets? Larger = bad
% compare is an array of scans that we compare to template scan

cumulative_offset = 0;

scans = [template compare];

for scanN=1:length(compare)
    
    % Importing data
    scan = compare(scanN);
    
    % Padding (to make maps the same size)
    if size(template.data,3) ~= size(scan.data,3)
        while size(template.data,3) > size(scan.data,3)
            scan.data = padarray(scan.data, [0 0 1 0], 'pre');
            fprintf('Added one row of padding to scan %d \n ', compare(scanN).data_meta.scan_number);
        end
        while size(template.data,3) < size(scan.data,3)
            template.data = padarray(template.data, [0 0 1 0], 'pre');
            fprintf('Added one row of padding to scan %d \n', template.data_meta.scan_number);
            for i=scanN-1:-1:1 %pads already aligned scans
                compare(i).data = padarray(compare(i).data, [0 0 1 0], 'pre');
            end
        end
    end
    if size(template.data,4) ~= size(scan.data,4)
        while size(template.data,4) > size(scan.data,4)
            scan.data = padarray(scan.data, [0 0 0 1], 'pre');
            fprintf('Added one column of padding to scan %d \n', scanN);
        end
        while size(template.data,4) < size(scan.data,4)
            template.data = padarray(template.data, [0 0 0 1], 'pre');
            fprintf('Added one column of padding to scan %d \n', template.data_meta.scan_number);
            for i=scanN-1:-1:1 
                compare(i).data = padarray(compare(i).data, [0 0 0 1], 'pre');
            end
        end
    end
    
    template.integrate([1 2]);
    scan.integrate([1 2]); 
    
    template_data = template.data_integral;
    compare_data = scan.data_integral;
    
%     % Calculation of offset
    corr = normxcorr2(template_data,compare_data);
    
    [max_value, max_index] = max(abs(corr(:)));
    [xpeak, ypeak] = ind2sub(size(corr),max_index(1));
    
    offset = [(xpeak-size(compare_data,1)) 
                   (ypeak-size(compare_data,2))];
    offset_x(scanN) = -offset(1);
    offset_y(scanN) = -offset(2);
    cumulative_offset = cumulative_offset + abs(offset_x(scanN)) + abs(offset_y(scanN));
    
    % Truncating data outside the window of overlay
                
    % x
    if offset_x(scanN) ~= 0
        
        cut_pixels_x = zeros(size(template_data,2),size(template_data,1),2);
        
        % Cut template
        if offset_x(scanN) > 0
            for i=1:offset_x(scanN)
                cut_pixels_x(:,i,1) = -1; 
            end
        end
        
        if offset_x(scanN) < 0
            for i=size(template_data,1):-1:size(template_data,1)+offset_x(scanN)+1
                cut_pixels_x(:,i,1) = -1; 
            end
        end
        
        template.cut(cut_pixels_x(:,:,1),'x');
        for i=scanN-1:-1:1 % cuts previously aligned scans too
            compare(i).cut(cut_pixels_x(:,:,1),'x');
        end
        
        % Cut scan
        offset_x(scanN) = -offset_x(scanN);
        
        if offset_x(scanN) > 0
            for i=1:offset_x(scanN)
                cut_pixels_x(:,i,2) = -1; 
            end
        end
        
        if offset_x(scanN) < 0
            for i=size(compare_data,1):-1:size(compare_data,1)+offset_x(scanN)+1
                cut_pixels_x(:,i,2) = -1; 
            end
        end
        
        scan.cut(cut_pixels_x(:,:,2),'x');
    end
    
    % y
    if offset_y(scanN) ~= 0
        % Cut template
        cut_pixels_y = zeros(size(template_data,2),size(template_data,1),2);
        
        if offset_y(scanN) > 0
            for i=1:offset_y(scanN)
                cut_pixels_y(i,:,1) = -1; 
            end
        end
        
        if offset_y(scanN) < 0
            for i=size(template_data,2):-1:size(template_data,2)+offset_y(scanN)+1
                cut_pixels_y(i,:,1) = -1; 
            end
        end
        
        template.cut(cut_pixels_y(:,:,1),'y');
        for i=scanN-1:-1:1 % cuts previously aligned scans too
            compare(i).cut(cut_pixels_y(:,:,1),'y');
        end
        
        % Cut scan
        offset_y(scanN) = -offset_y(scanN);
        
        if offset_y(scanN) > 0
            for i=1:offset_y(scanN)
                cut_pixels_y(i,:,2) = -1; 
            end
        end
        
        if offset_y(scanN) < 0
            for i=size(compare_data,2):-1:size(compare_data,2)+offset_y(scanN)+1
                cut_pixels_y(i,:,2) = -1; 
            end
        end
        
        scan.cut(cut_pixels_y(:,:,2),'y');
    end
end

%Cut the padding
for scanN=1:length(scans)
    
    scan = scans(scanN);
    
    if scanN==1
        scan.integrate([1 2]); 
        cut_pixels = zeros(size(scan.data_integral,1),size(scan.data_integral,2));
    
       
        [I,J] = ind2sub(size(scan.data_integral),find(scan.data_integral==0));
        
        for row=1:length(I)
            for col=1:length(J)
                cut_pixels(I(row),J(col)) = -1; 
            end
        end
        
        cut_pixels = transpose(cut_pixels);
        
        col_cut = all(cut_pixels==-1, 1);
        cut_pixels_x = zeros(size(cut_pixels));
        
        for i=1:length(col_cut)
            if col_cut(i) == 1
                cut_pixels_x(:,i) = -1;
            end
        end
        
        row_cut = all(cut_pixels==-1, 2);
        cut_pixels_y = zeros(size(cut_pixels));
        
        for i=1:length(row_cut)
            if row_cut(i) == 1
                cut_pixels_y(i,:) = -1;
            end
        end
    end
    
    scan.cut(cut_pixels_x,'x'); 
    scan.cut(cut_pixels_y,'y');
    
%     log
    aligned_with = '';
    for j=1:length(scans)
        if scanN ~= j
            aligned_with = [aligned_with ' ' num2str(scans(j).data_meta.scan_number)];
        end
    end
    scan.log = [scan.log 'Aligned with scans' aligned_with ' '];
    scan.reset;
    
end

disp('Scans were aligned. ');
