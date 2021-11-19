% ********************************************************************** %
% Preprocessing Script for DEED Lab Resting State EEG Data [Script 1]
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

% Path of folder with raw data
% Data is in .mff (EGI/Magstim) format
proj.data_location = 'INSERT_PATH_WITH_FOLDER';

% Get mff file names
proj.mff_filenames = dir(fullfile(proj.data_location, '*.mff'));
proj.mff_filenames = { proj.mff_filenames(:).name };

% Location for a file to hold error messages for subjects whose processing fails
proj.error_file = 'INSERT_PATH_WITH_FOLDER\errors.txt';

% ********************************************************************** %

%% Loop over subjects and run rest_process_single_subject.m

for i = 1:length(proj.mff_filenames)
    proj.currentSub = i;
    proj.currentId = proj.mff_filenames{i};
    
    % Subject ID will be filename up to first space, or up to first '.'
    space_ind = strfind(proj.currentId, ' ');
    if ~isempty(space_ind)
        proj.currentId = proj.currentId(1:(space_ind(1)-1)); 
    else
        mff_ind = strfind(proj.currentId, '.mff');
        proj.currentId = proj.currentId(1:(mff_ind(1)-1));
    end
    
    try
        if i == 1
            summary_info = rest_process_single_subject;
            summary_tab = struct2table(summary_info);
        else
            summary_info = rest_process_single_subject;
            summary_row = struct2table(summary_info); % 1-row table
            summary_tab = vertcat(summary_tab, summary_row); % Append new row to table
        end
        
    catch me
       fid = fopen(proj.error_file, 'a');
       % At (date) x at time y subject z had error q
       fprintf(fid, 'At %s subject %s had error %s\r\n', ...
           datestr(now), proj.currentId, me.message);
       fprintf(fid, '\tin %s at line %s \r\n', me.stack(end-1).file, ...
           num2str(me.stack(end-1).line)); 
       fclose(fid);
    end
    
end

%% Write summary info to spreadsheet

proj.output_location = 'INSERT_PATH_WITH_FOLDER';
writetable(summary_tab, [proj.output_location filesep 'preprocessing_log.csv']);

% ****************************** THE END ******************************* %