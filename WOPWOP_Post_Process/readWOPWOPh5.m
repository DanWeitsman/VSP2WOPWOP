function [data] = readWOPWOPh5(filename)

%This function converts the contents of an HDF5 file to a matlab struct.
%The output struct will contain 3 fields; functional_values, geometry_values, and names.
%'function_values' contains the predicted pressure time series for each
%observer and is sized as (time,thickness,loading,total) x (time series) x(# of observers).
% 'geometric_values' provides the coordinates of each observer but is sized just as 'functional_values',
%meaning the second dimension of this array is redundant. The 'names' field just provide the names of each 
%of the data sets found in the first dimension of the 'function_values' field. 

info = h5info(filename);
dset_name = {info.Datasets.Name};

for i = 1:length(dset_name)
    data.(dset_name{i}) = h5read(filename,strcat('/',dset_name{i}));
end
end

