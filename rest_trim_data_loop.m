% ********************************************************************** %
% Preprocessing Script for DEED Lab Resting State EEG Data [Script 4]
% Authors: Armen Bagdasarov (Graduate Student) & 
%          Kenneth Roberts (Research Associate)
% Institution: Duke University
% Year Created / Last Modified: 2021 / 2021
% ********************************************************************** %

%% Prepare workspace for preprocessing

% Clear workspace and command window
clear;
clc;

% Start EEGLAB 
% Startup file in MATLAB folder should have already added it to the path
eeglab;

% Declare variables as global
% Global variables are those that you can access in other functions
global proj;

% Path of folder with data in .set format
proj.data_location = 'INSERT_PATH_WITH_FOLDER';

% Get set file names
proj.set_filenames = dir(fullfile(proj.data_location, '*.set'));
proj.set_filenames = { proj.set_filenames(:).name };

%% Loop over subjects and run rest_trim_data_single_subject.m

for i = 1:length(proj.set_filenames)
    proj.currentSub = i;
    proj.currentId = proj.set_filenames{i};
    
    % Subject ID will be filename up to first underscore, or up to first '.'
    space_ind = strfind(proj.currentId, '_');
    if ~isempty(space_ind)
        proj.currentId = proj.currentId(1:(space_ind(1)-1)); 
    else
        set_ind = strfind(proj.currentId, '.set');
        proj.currentId = proj.currentId(1:(set_ind(1)-1));
    end
    
    if i == 1
        summary_info = rest_trim_data_single_subject;
        summary_tab = struct2table(summary_info);
    else
        summary_info = rest_trim_data_single_subject;
        summary_row = struct2table(summary_info); % 1-row table
        summary_tab = vertcat(summary_tab, summary_row); % Append new row to table
    end
    
end

% ****************************** THE END ******************************* %