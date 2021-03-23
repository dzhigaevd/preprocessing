function [binned] = bin(image, binning_size)
    convoluted = conv2(image, ones(binning_size));
    convoluted_size = size(convoluted);
    binned = convoluted(binning_size:binning_size:convoluted_size(1), binning_size:binning_size:convoluted_size(2));
end