% ********************************************************************** %
% Preprocessing Script for DEED Lab Resting State EEG Data [Script 3]
% Authors: Armen Bagdasarov (Graduate Student) & 
%          Kenneth Roberts (Research Associate)
% Institution: Duke University
% Year Created / Last Modified: 2021 / 2021
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

function [EEG, info] = create_eyes_open_closed_resting_events(EEG)

% From resting files with 8 'bgin' events, this will:
%   1) Relabel those events with rs_open and rs_closed
%   2) Generate 'rs_end' events to mark 60s after beginning of resting state
%   3) Do checks and warn if any epochs overlap
%   4) Do checks and warn if any epoch is too short (goes off end of file)

% For debug/development purposes, if first argument in is string 'debug',
% then do it with this hard-coded file
persistent old_EEG;
if ~isstruct(EEG) && strcmp(EEG, 'debug')
    if isempty(old_EEG)
        test_filename = 'PATH_WITH_FILE';
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        eeglab('redraw');
        EEG = pop_mffimport({test_filename},{'code'});
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
        old_EEG = EEG;
    else
        EEG = old_EEG;
    end
end

% Read data out of event
evt_codes = { EEG.event(:).code };
evt_types = { EEG.event(:).type };
evt_lats = [ EEG.event(:).latency ];

% Find 'bgin' events
rs_evt_ind = find(strcmp('bgin', evt_codes));

% Label 'open' and 'closed'
[ EEG.event(rs_evt_ind(1:2:end)).code ] = deal('rs_open');
[ EEG.event(rs_evt_ind(1:2:end)).type ] = deal('rs_open');
[ EEG.event(rs_evt_ind(2:2:end)).code ] = deal('rs_closed');
[ EEG.event(rs_evt_ind(2:2:end)).type ] = deal('rs_closed');

% Add new events that mark the end of each period
nevt = length(EEG.event);
evt_template = EEG.event(rs_evt_ind(1));   
    % Make a "template" event from first rs event
evt_newlats = [EEG.event(rs_evt_ind).latency] + 60*EEG.srate;

% Ensure no event goes off end of the file, move events to within 10 pnts
% of end of file (leaving a little buffer for resampling/filtering)
evt_newlats_from_end = EEG.pnts - evt_newlats;
evt_newlats = min(evt_newlats, EEG.pnts-10);

[~, si] = sort([evt_lats, evt_newlats]);    
    % Need to sort in case a block is shorter than 60s
EEG.event(si <= nevt) = EEG.event(:);       
    % Spread out old events

% Fill in values of newm 'rs_end' events
[EEG.event(si > nevt)] = deal(EEG.event(rs_evt_ind(1)));    
    % Generate new events by copying first rs 'bgin'
[EEG.event(si > nevt).code] = deal('rs_end');               
    % Fix events by filling in important fields
[EEG.event(si > nevt).type] = deal('rs_end');
[EEG.event(si > nevt).begintime] = deal('');

evt_newlats = mat2cell(evt_newlats, [1], ones(1, length(evt_newlats)));
[EEG.event(si > nevt).latency] = evt_newlats{:};

% See if blocks overlap
begin_ind = ismember({EEG.event(:).code}, {'rs_open', 'rs_closed'});
end_ind = ismember({EEG.event(:).code}, {'rs_end'});

blocklen = diff([EEG.event(begin_ind).latency EEG.pnts]/EEG.srate);
short_block = find(blocklen < 60);
if ~isempty(short_block)
    warning('Blocks %s too short, length of %s seconds respectively.\n', ...
        mat2str(short_block), mat2str(blocklen(short_block)));
    info.block_overlap = true;
    info.block_int = blocklen;
else
    info.block_overlap = false;
    info.block_int = blocklen;
end

% See if any block has been ended early because it exceeded length of file 
% or came within 10 points
if any(evt_newlats_from_end < 10)
    warning('Final %d blocks were truncated by early end of file.', ...
        sum(evt_newlats_from_end <0));
    info.block_truncate = true;
else
    info.block_truncate = false;
end

info.blocklen = ([EEG.event(end_ind).latency] - [EEG.event(begin_ind).latency])/EEG.srate;

% ****************************** THE END ******************************* %