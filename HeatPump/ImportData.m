%% Import data from text file.
% Script for importing data from the following text file:
%
%    /home/thomas/ownCloud2/MES3/Semester project/Code/heat_data2019_07_01to2020_01_31.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2020/04/07 12:41:21

%% Initialize variables.
filename = '/home/thomas/ownCloud2/MES3/Semester project/Code/heat_data2019_07_01to2020_01_31.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: double (%f)
%	column2: datetimes (%{yyy-MM-dd HH:mm:ss}D)
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%{yyy-MM-dd HH:mm:ss}D%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
heatdata = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','time','heat_consumption'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% heatdata20190701to1.time=datenum(heatdata20190701to1.time);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;