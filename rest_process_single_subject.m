% ********************************************************************** %
% Preprocessing Script for DEED Lab Resting State EEG Data [Script 2]
% Authors: Armen Bagdasarov (Graduate Student) & 
%          Kenneth Roberts (Research Associate)
% Institution: Duke University
% Year Created / Last Modified: 2021 / 2021
% ********************************************************************** %

%% Run loop_over_subjects script to run this function

function summary_info = rest_process_single_subject(varargin)
%% Declare variables as global
% Global variables are those that you can access in other functions
global proj

% ********************************************************************** %

%% Reset the random number generator to try to make the results replicable 
% This produces somewhat consistent results for clean_rawdata functions 
% and ICA which can vary by run
rng('default');

% ********************************************************************** %

%% Import data
% Raw data is in .mff format (Magstim/EGI)
mff_filename = fullfile(proj.data_location, ...
    proj.mff_filenames{proj.currentSub}); % Get file name
EEG = pop_mffimport({mff_filename},{'code'}); % Import .mff
summary_info.currentId = {proj.currentId}; % Save subject ID in summary info

% ********************************************************************** %

%% Remove outer ring of electrodes
% The outer ring is often noisy in high-density nets, 
% so we remove these channels

% Outer layer of channels to be removed 
% 25 of them, including the online reference E129/Cz which is flat
outer_chans = {'E129' 'E17' 'E38' 'E43' 'E44' 'E48' 'E49' 'E113' 'E114' ...
    'E119' 'E120' 'E121' 'E125' 'E126' 'E127' 'E128' 'E56' 'E63' 'E68' ...
    'E73' 'E81' 'E88' 'E94' 'E99' 'E107'};

% Remove outer_chans
EEG = pop_select(EEG, 'nochannel', outer_chans);

% Save variable with reduced (104) channel locations 
% This will be needed later when interpolating bad channels
reduced_chan_locs = EEG.chanlocs;

% ********************************************************************** %

%% Downsample from 1000 to 250 Hz 
EEG = pop_resample(EEG, 250);

% Note: This step must be done before inserting event markers
% Otherwise, it produced an error (not sure why, does not matter)

% ********************************************************************** %

%% Insert eyes- open and closed markers
% The current data does not contain markers for when the condition is 
% eyes- open vs. closed; only markers for when a condition begins
% Need to insert markers
% 8 runs, each 60 seconds or 1 minute long in the following order: 
    % 1. Open 
    % 2. Closed 
    % 3. Open
    % 4. Closed
    % 5. Open
    % 6. Closed
    % 7. Open
    % 8. Closed

% This following is a function from another script [Script 3]
[EEG, info] = create_eyes_open_closed_resting_events(EEG);

% Save whether any blocks overlap (binary 0/1 = no/yes) in summary info
% This would be bad
summary_info.block_overlap = info.block_overlap;

% Save whether any block ended early (binary 0/1 = no/yes) in summary info
% This would be bad
summary_info.block_truncate = info.block_truncate;

% Keep only the eyes-closed data (in other words remove the eyes-open data)
% After this step, data will be 240 seconds
EEG = pop_rmdat(EEG, {'rs_closed'}, [0 60] ,0);

% Save length of data as an additional check
summary_info.data_length_check_240 = EEG.xmax;

% ********************************************************************** %

%% Filter & Remove Line Noise

% Low-pass filter at 40 Hz and high-pass filter at 1 Hz
EEG = pop_eegfiltnew(EEG, 'hicutoff', 40); 
    % Because low-pass filter is applied at 40 Hz (needed for microstates), 
    % don't have to worry too much about line noise at 60 Hz
    % However, there will still be some line noise
EEG = pop_eegfiltnew(EEG, 'locutoff', 1); 
    % 1 Hz is necessary for ASR, ICA, and microstates
    
% CleanLine to remove 60 Hz line noise
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan ,'computepower',1,...
    'linefreqs',60,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,...
    'plotfigures',0,'scanforlines',1,'sigtype','Channels','taperbandwidth',2,...
    'tau',50,'verb',1,'winsize',4,'winstep',1);
% All are default parameters, except tau changed from 100 to 50, 
% which is recomended by the authors for continuous unepoched data

% ********************************************************************** %

%% Reject bad channels

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off',...
    'WindowCriterion','off','BurstRejection','off','Distance','Euclidian',...
    'MaxMem', 48); 
        % All default settings
        % Only rejecting bad channels, everything else (ASR) turned off
        % MaxMem = 48gb for reproducibility (will vary based on RAM)

% Save which channels were bad in summary info
if isempty(EEG.chaninfo.removedchans) % If no bad chans...
    summary_info.bad_chans = {[]}; % ...then leave blank
else
    bad_chans = {EEG.chaninfo.removedchans(:).labels};
    summary_info.bad_chans = {strjoin(bad_chans)};

    % Plot bad channels to identify whether there are clusters of bad channels
    bad_chan_ind = find(ismember({reduced_chan_locs(:).labels}, bad_chans));
    figure; topoplot(bad_chan_ind, reduced_chan_locs, 'style', 'blank', ...
        'emarker', {'.','k',[],10}, 'electrodes', 'ptslabels');

    % Save bad channels plot
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
        'Inches', 'PaperSize', [10, 10])
    bad_chan_plot_path = 'INSERT_PATH_WITH_FOLDER'; 
    bad_chan_plot_name = [proj.currentId '_rest_bad_chans_plot'];
    saveas(gca, fullfile(bad_chan_plot_path, bad_chan_plot_name), 'png');
    close(gcf);
end

% Save number of bad channels in summary info
summary_info.n_bad_chans = length(EEG.chaninfo.removedchans);

%% Remove large artifacts

% Artifact Subspace Reconstruction (ASR) + 
% additional removal of bad data periods

% First, save data before ASR
% ICA will be run later on the data post-ASR
% But ICA fields will be applied to the to the pre-ASR data
EEG_no_rej = EEG;

% ASR
% All default settings
% Most importantly the burst criterion is set conservatively to 20 
% and burst rejection is set to on (meaning remove instead of fix bad data)
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off',...
    'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
    'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
    'WindowCriterionTolerances',[-Inf 7], 'MaxMem', 48); 
        % MaxMem set to 48gb for reproducibility 

% Save how many seconds of data is left after ASR in summary info
% This will be important later for excluding participants 
% For example, if ICA was run only on 30 seconds of data because ASR cut
% out the rest, should get rid of the file (i.e., their
% data was probably very noisy)... not enough data for ICA to be reliable
summary_info.post_ASR_data_length = EEG.xmax;

% ********************************************************************** %

%% ICA
% Extended infomax ICA with PCA dimension reduction
% PCA dimension reduction is necessary because of the large number of 
% channels and relatively short amount of data
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', 30);

% Save ICA plot
ica_plot_path = 'INSERT_PATH_WITH_FOLDER'; 
ica_plot_name = [proj.currentId 'rest_ica_plot'];
pop_topoplot(EEG, 0, [1:30], 'Independent Components', 0, 'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
    'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_plot_path, ica_plot_name), 'png');
close(gcf);

% ********************************************************************** %

%% Select independent components (ICs) related to eye or muscle only

% Automatic classification with ICLabel
EEG = pop_iclabel(EEG, 'default');

% Flag components with >= 70% of being eye or muscle
EEG = pop_icflag(EEG, [NaN NaN; 0.7 1; 0.7 1; NaN NaN; NaN NaN; ...
    NaN NaN; NaN NaN]);

% Select components with >= 70% of being eye or musscle
eye_prob = EEG.etc.ic_classification.ICLabel.classifications(:,3);
muscle_prob = EEG.etc.ic_classification.ICLabel.classifications(:,2);
eye_rej = find(eye_prob >= .70);
muscle_rej = find(muscle_prob >= .70);
eye_muscle_rej = [eye_rej; muscle_rej];
eye_muscle_rej = eye_muscle_rej';

% Save retained variance post-ICA in summary info
[projected, pvar] = compvar(EEG.data, {EEG.icasphere, EEG.icaweights}, EEG.icawinv, eye_muscle_rej);
summary_info.var_retained = 100-pvar;

% Plot only the removed components
ica_rej_plot_path = 'INSERT_PATH_WITH_FOLDER'; 
ica_rej_plot_name = [proj.currentId '_rest_removed_ics'];
figure % This line is necessary for if there is only 1 component to plot
pop_topoplot(EEG, 0, eye_muscle_rej, 'Independent Components', 0, ...
    'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', ...
    'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_rej_plot_path, ica_rej_plot_name), 'png');
close(gcf); close(gcf); % Need to close twice if there is only 1 component to plot

% ********************************************************************** %
%% Copy EEG ICA fields to EEG_no_rej and remove ICs with >= 70% of being eye or muscle
% Basically, back-projecting the ICA information from the ASR-reduced data 
% to the full pre-ASR data

EEG_no_rej.icawinv = EEG.icawinv;
EEG_no_rej.icasphere = EEG.icasphere;
EEG_no_rej.icaweights = EEG.icaweights;
EEG_no_rej.icachansind = EEG.icachansind;

EEG = EEG_no_rej; % Set EEG to the one with full data length, pre-ASR

% Remove components with >= 70% of being eye or muscle
EEG = pop_subcomp(EEG, eye_muscle_rej , 0);

if isempty(eye_muscle_rej) % If no ICs removed...
    summary_info.ics_removed = {[]}; % ...then leave blank
else
    % Save which components were removed in summary info
    summary_info.ics_removed = {num2str(eye_muscle_rej)};
end

% Save the number of components removed in summary info
summary_info.n_ics_removed = length(eye_muscle_rej);

% ********************************************************************** %

%% Additional artifact rejection
% At this point, the file is still the original length (240 seconds)
% But it likely contians artifact 
% Epoch the data and use the TBT plug-in to reject epochs
% This plug-in also does epoch-by-epoch channel interpolation, which is nice

% First epoch the data into 1 second segments
% There are some discontinuities so will not get exactly 240 epochs
% Instead 237 epochs; lose 3 epochs / 3 seconds of data here (that's okay)
EEG = eeg_regepochs(EEG, 'recurrence', 1, 'rmbase', NaN);

% 1. Abnormal values
% Simple voltage thresholding of -100/+100 uV
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -100 , 100 ,0 , 0.996, 1, 0);

% 2. Improbable data
% Based on joint probability, SD = 3 for both local and global thresholds 
EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 3, 3, 1, 0, 0);

% Reject based on the epochs selected above
EEG = pop_TBT(EEG, EEG.reject.rejthreshE | EEG.reject.rejjpE, 10, 1, 0);
    % Do epoch interpolation on both types of artifact rejection at once
    % Criteria must be met in at least 10 channels for the epoch to be rejected
    % Don't want to remove channels, so criteria must be met in all 
    % channels for the channel to be removed, which is unlikely to happen

% Save how many epochs are left 
summary_info.n_epochs = EEG.trials;

% ********************************************************************** %

%% Interpolate removed bad channels
EEG = pop_interp(EEG, reduced_chan_locs, 'spherical');

% ********************************************************************** %

%% Re-reference to the average
EEG = pop_reref(EEG, []);

% ********************************************************************** %

%% Plot Channel Spectra

% Plot channel spectra 
spectra_time = EEG.xmax * 1000;
figure; pop_spectopo(EEG, 1, [0  spectra_time], 'EEG' , 'freqrange',...
    [2 80],'electrodes','off');

% Save channel spectra plot
spectra_plot_path = 'INSERT_PATH_WITH_FOLDER'; 
spectra_plot_name = [proj.currentId '_rest_spectra_plot'];
saveas(gca, fullfile(spectra_plot_path, spectra_plot_name), 'png');
close(gcf);

% ********************************************************************** %

%% Save final preprocessed files
% in .set format
set_path = 'INSERT_PATH_WITH_FOLDER'; 
set_name = [proj.currentId '_rest_complete'];
pop_saveset(EEG, fullfile(set_path, set_name));

% ****************************** THE END ******************************* %
