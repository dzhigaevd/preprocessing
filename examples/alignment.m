%% Align drift
clc;
clear all;

scan_numbers = 185:215;

for i=1:length(scan_numbers)
    path = sprintf('/home/jlastow/Documents/data/analysis/sample_GaSb_W/Sample1_%d/Sample1_%d_merlin.mat',scan_numbers(i),scan_numbers(i));
    load(path);
    scans(i) = obj;
    disp('Loaded scan.')
    scans(i).drift_align('x');
    scans(i).save_data;
    i
end

%% Align angles

%clear all;
clc;

scan_numbers = 185:215

% Loading
for i=1:length(scan_numbers)
    path = sprintf('/home/jlastow/Documents/data/analysis/sample_GaSb_W/Sample1_%d/Sample1_%d_merlin.mat',scan_numbers(i),scan_numbers(i));
    load(path);
    if i == 1
        template(i) = obj
    else
        compare(i-1) = obj
    end
    disp('Loaded scan.')
    i
end

% Alignment
cumulative_offset = align_scans(template, compare)

% Saving data
for i=1:length(scan_numbers)
    if i == 1
        template(i).save_data;
    else
        compare(i-1).save_data;
    end
    i
end

%% Construct 5D matrix

clc;
clear all;

scan_numbers = 185:215;

for i=1:length(scan_numbers)
    path = sprintf('/home/jlastow/Documents/data/analysis/sample_GaSb_W/Sample1_%d/Sample1_%d_merlin.mat',scan_numbers(i),scan_numbers(i));
    load(path);
    %scans(i) = obj;
    fprintf('Loaded scan %d/%d', i, length(scan_numbers))
    data(:,:,:,:,i) = obj.data;
end

save('data185to215', 'data', '-v7.3')

%% Rocking curve

scan_numbers = 184:215;
integrated = []

for i=1:length(scan_numbers)
    if i == 1
        template(i).integrate([1 2 3 4]);
        integrated(i) = template(i).data_integral;
    else
         compare(i-1).integrate([1 2 3 4]);
        integrated(i) = compare(i-1).data_integral;
    end
end

plot(integrated)

%% Flux

scan_numbers = 202:208;

for i=1:length(scan_numbers)
    if i == 1
        flux_i = openh5attribute(template(i).data_meta.master_file_nanomax, sprintf('/entry%d/measurement/Ni6602_buff',template(i).data_meta.scan_number));
        nbr = template(i).data_meta.scan_number;
    else
        flux_i = openh5attribute(compare(i-1).data_meta.master_file_nanomax, sprintf('/entry%d/measurement/Ni6602_buff',compare(i-1).data_meta.scan_number));
        nbr = compare(i-1).data_meta.scan_number;
    end
    flux(:,:,i) = flux_i;
    figure();
    imagesc(flux_i')
    colorbar;
    title(nbr)
    flux_val(i) = squeeze(sum(flux_i,[1 2])); %integrated value
end

figure();
plot(flux_val)

