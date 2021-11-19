% ********************************************************************** %
% Preprocessing Script for DEED Lab Resting State EEG Data [Script 5]
% Authors: Armen Bagdasarov (Graduate Student) & 
%          Kenneth Roberts (Research Associate)
% Institution: Duke University
% Year Created / Last Modified: 2021 / 2021
% ********************************************************************** %

%% Run loop_over_subjects script to run this function

function summary_info = rest_trim_data_single_subject(varargin)
%% Declare variables as global
% Global variables are those that you can access in other functions
global proj;

% ********************************************************************** %

%% Import data
set_filename = proj.set_filenames{proj.currentSub};
EEG = pop_loadset('filename', {set_filename}, 'filepath',...
    'INSERT_PATH_WITH_FOLDER');
summary_info.currentId = {proj.currentId};

% ********************************************************************** %

%% Trim data to 145 seconds / 145 epochs
EEG = pop_select( EEG, 'trial',[1:145]);

% ********************************************************************** %

%% Save final preprocessed files in .set and .bva formats

% .set format
set_path = 'INSERT_PATH_WITH_FOLDER'; 
set_name = [proj.currentId '_' 'rest' '_' 'processed' '_' '145s'];
pop_saveset(EEG, fullfile(set_path, set_name));

% .bva format (for CARTOOL import)
bva_path = 'INSERT_PATH_WITH_FOLDER'; 
bva_name = [proj.currentId '_' 'rest' '_' 'processed' '_' '145s'];
pop_writebva(EEG, fullfile(bva_path, bva_name));

% ****************************** THE END ******************************* %