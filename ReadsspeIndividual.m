function [dat] = ReadsspeIndividual(file)
%
dirFunctionLib = 'N:\4all\public\mpsd_dmiller\Stuart\MATLAB\data_processing_libs\Data_Processing2018';
addpath(dirFunctionLib);
dataIOAssembly = NET.addAssembly([dirFunctionLib,'\DataIOLib.dll']);

% Individual backup frames are saved as 'uint16' format
% fid=fopen(file); 
% Im1=fread(fid,'int32'); % reads the entire file, including metadata
% 
% px=Im1(2);
% py=Im1(3);
% 
% Im=Im1(4:end);
% Z=reshape(Im,px,py);
% Z=reshape(Im1(2:end-4),1280,1080);

dat=double(DataIOLibrary.DataIO.ReadSpe(file));

% fclose(fid);
end