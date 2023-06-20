%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script for data pre-processing of the infant audiovisual habituation eyetracking experiment (SGU Task 2 - 2023).
% 
% Data collected March 2023.
% 
% This script was originally coded for the infant audiovisual task (2022) by Ana Maria Portugal then adapted for
% this project by Giorgia Bussu.
%
% The script uses the hdf5 file and the trial information csv file.
%
% The task script was implemented by Giorgia Bussu in Python (Psychopy).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all


%% Set general variables

% for interpolation, both pupil and gaze, the maximum gap allowed
interp_maxGap = 150; % in ms, leave empty if no interpolation is needed

% path for stimuli that define the stimulus condition in the csv trial file
% (unisensory vs. multisensory event)
multi = 'stimuli/habituation_videos/crash/';
uni = 'stimuli/habituation_videos/no_crash/';


% Screen limits in the gaze data (distance to screen assumed to be 65 cm)
% Monitor size is 52.69x29.64. Monitor right-center edge has 
% coordinates (22.06, 0), the  center-top edge has coordinates (0, 12.9). 
% The left-bottom corner of the monitor has coordinates (-22.06, -12.9).
width_half = 22.06; % X from -22.06 to 22.06
height_half = 12.9; % Y from -12.9 to 12.9

% Define AOI size - during task presentation the stimuli was within a 
% 28x21 visual degrees size
stim_box_grab = 11;
stim_box = 28; % AOI is right 50% screen

% do we want plots and data structs to be saved?
plot_figure = 1; % 1=saves heatmap and a figure with 4 different plots 
save_data = 1; % 1=saves wide and long csv files

% initiate a matrix to store pupil samples across trials and individuals
t_for_B = 1; %this is needed to collect all pupil samples in the baseline 
t_for_R = 1; %this is needed to collect all pupil samples in response 


%% paths for data

resultsPath = 'C:\Users\giobu365\Documents\pupil_habituation\ET_DATA';
resultsPath_calib = 'C:\Users\giobu365\Documents\pupil_habituation\ET_calibration';

cd(resultsPath)
mkdir('..\results_habituation_shortwindow2\') % make directory for plots and results
files = dir('**\*.csv'); % scan folder for all trial info csv files, 
% should be one csv file per participant

for f = 1: length(files)
    
    if f > length(files)
        % this is needed because we'll remove some participants below
        continue
    end
    
    [~,name,ext] = fileparts(files(f).name);
    subject = name(1:5); % based on experiment subject codes 
   
    s_time = name(end-3:end);
    id_name = [subject, '_', s_time]; % time was saved because we had duplicates
    
%     % uncomment this to inspect particular subjects
%         if subject == "P17"
%         else
%             continue
%         end
        
    
    %% skip/ remove duplicates - this needs adaptation in case of other projects
    % when specific subjects and times appear we removed their paths from 
    % the files list and reassigned the subject variables
    
    if subject == "SGU22" && s_time == "1110"
        warning('skip')
        files(f)=[];
        if f > length(files) % if we are at the end of the files list we
            % move on in the loop / end loop
            continue
        else % if not we re assign the subject name and paths to the loop index            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end
    end
    if subject == "SGU22" && s_time == "1113"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end
    end
    if subject == "SGU27" && s_time == "1531"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end
    end
    if subject == "SGU27" && s_time == "1535"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end
    end
    if subject == "SGU27" && s_time == "1552"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end
    end
    if subject == "SGU29" && s_time == "1203"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end 
    end
    if subject == "SGU29" && s_time == "1210"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:5);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time];
        end 
         
    end
    
    % trial_data_all_wide collects wide format data, one row per subject
    trial_data_all_wide.subject(f,1) = convertCharsToStrings(subject);
    trial_data_all_wide.id_name(f,1) = convertCharsToStrings(id_name);
            
    
    %% get invidiual file paths
    
    % find the hdf5 file - this file has the same filename as the csv file
    % but ends in '_hdf5.hdf5'
    Filename = dir(['**/',name,'_hdf5.hdf5']); 
    if isempty(Filename)
        % if no file was found then move on to the next subject
        disp(warning('No hdf5 file found'))
        continue
    else
        Filename =  [Filename.folder,'\', Filename.name ];
    end
    
    % find the trial information csv file - this just needs to get the
    % pathname and filename already in the files struct
    Filename_csv = [files(f).folder,'\', files(f).name ];
 
    % find the calibration file - this searches based on the subject id 
    % it might give us multiple files or files that are not necessarily
    % corresponding to the correct session. This was fixed after manually
    Filename_calib = dir([resultsPath_calib, '\',subject,'*ease_et_calibration*', '_validation_data.csv' ]);
    
    if isempty(Filename_calib)
        % if no file was found, write NaN in the variables field
        disp(warning('No Calibration file found'))
        trial_data_all_wide.calib_time(f,1) = "No file";
        
        trial_data_all_wide.mean_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.std_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.proportion_of_time_gaze_on_screen(f,1) = NaN;

    elseif length(Filename_calib) > 1
        % if more than on file was found, write NaN in the variables field
        disp(warning('More than one calibration file found'))
        trial_data_all_wide.calib_time(f,1) = "More than one file, fix manually";
        
        trial_data_all_wide.mean_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.std_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.proportion_of_time_gaze_on_screen(f,1) = NaN;
        
    else   
        % if we find only one file, we will get the time on the file so we
        % can check if the corresponds to the dataset (i.e. should be some
        % minutes earlier than the eye tracking session
        trial_data_all_wide.calib_time(f,1) = convertCharsToStrings(Filename_calib.name(end-33:end-20));
        
        Filename_calib =  [Filename_calib.folder,'\', Filename_calib.name ];
        calib_info = ease_2_getcalibinfo(Filename_calib); 
        
        trial_data_all_wide.mean_distance_gaze_to_target(f,1) = calib_info.mean_distance_gaze_to_target;
        trial_data_all_wide.std_distance_gaze_to_target(f,1) = calib_info.std_distance_gaze_to_target;
        trial_data_all_wide.proportion_of_time_gaze_on_screen(f,1) = calib_info.proportion_of_time_gaze_on_screen;        
    end
    
    
    %% load gaze and event buffer from hdf5 file
    gaze_buffer  = h5read(Filename,'/data_collection/events/eyetracker/BinocularEyeSampleEvent');
    event_buffer = h5read(Filename, '/data_collection/events/experiment/MessageEvent');
   
    test  = h5read(Filename,'/data_collection/session_meta_data'); % this only has the date of the session
    experiment_date = string(test.code(3:18)');
    
    %test  = h5read(Filename,'/data_collection/events/experiment/LogEvent'); % this only has system information
    
    % convert char event information to strings for easy reading and access
    for i=1: size(event_buffer.text,2)
        idx= find(isletter(event_buffer.text(:, i)), 1, "last" );
        event_buffer.event(i,1) = convertCharsToStrings(event_buffer.text(1:idx, i)');
    end
        
    
    %% load trial information from the csv file
        
    trial_info = habituation_popout_gettrialinfo(Filename_csv);
    
    
    % remove rows that do not have task trial information / block breaks
    toRemove = find(isnan(trial_info.trial_habituation_global_start_time));
    trial_info(toRemove,:) = [];


    %% loop across all trials
    
    % choose which trials to plot if plot_figure = 1 above 
    trials_to_plot = randi([1 length(trial_info.trial_loop_thisRepN) ],1,2); % this plots two random trials
    %trials_to_plot = [1, 2, 3]; % this plots specific trials
    %trials_to_plot = [1:length(trial_info.trial_loop_thisRepN)]; % this plots all trials
  
    
    for t= 1:length(trial_info.trial_habit_loop_thisRepN)
        
%         % uncomment this to inspect particular trials
%             if t ~= 45   
%                 continue
%             end
        
        
        % trial_data collects long format data, one row per trial per
        % subject
        % save trial information - id and trial nr
        trial_data.id_name(t,1) = convertCharsToStrings(id_name);
        trial_data.id_date(t,1) = experiment_date;
        trial_data.id_trial(t,1) = t;
        
        %% save trial conditions - type, volume, and order
        
        if strncmp(string(trial_info.habituation_video_filepath(t)), multi, length(multi))
            trial_data.StimulusType(t,1) = 1; % multisensory event / crash
        elseif strncmp(string(trial_info.habituation_video_filepath(t)), uni, length(uni))
            trial_data.StimulusType(t,1) = 2; % unisensory event / no_crash
        end
        
        % we get the trial order using the trial number from the trial file
        trial_data.Order(t,1) = trial_info.trial_habit_loop_thisRepN(t);
        
       
        %% create AOIs for right half of the screen + central fixation (attention grabber)

        habit_AOI = [0, -height_half, width_half, 2*height_half]; %[ 0, - stim_box/2, stim_box/2, 1]; 
        fixation_AOI = [ -stim_box_grab/2, - stim_box_grab/2, stim_box_grab/2, stim_box_grab/2]; % central stim
                       
        
        %% segment trial
        
%         % compare global start time on csv to start of trial
%         trial_global_start_time = trial_info.trial_global_start_time(t);
%         trial_start_ebidx_global = find(event_buffer.event == ['exp1 trial ', num2str(t), ' start' ] );
%         trial_start_ebtime_global = event_buffer.time(trial_start_ebidx_global);
%         trial_global_start_time - trial_start_ebtime_global

        % trial start
        % the start of the trial will be 500 ms before the gaze was
        % captured / the offset of the attention grabber 
        % trial_start_ebtime might precede trial start but we want the 
        % window where gaze was evaluated to be in the fixation AOI

        if t<7
            ttrial = t+7;
            trial_data.Block(t,1) = 1;
        elseif t>6 && t<13
            ttrial = t+14;
            trial_data.Block(t,1) = 2;
        elseif t>12 && t<19
            ttrial = t+21;
            trial_data.Block(t,1) = 3;
        elseif t>18 && t<25
            ttrial = t+28;
            trial_data.Block(t,1) = 4;
        elseif t>24 && t<31
            ttrial = t+35;
            trial_data.Block(t,1) = 5;
        elseif t>30 && t<37
            ttrial = t+42;
            trial_data.Block(t,1) = 6;
        end

        trial_start_ebidx = find(event_buffer.event == ['exp2 habituation trial ', num2str(ttrial), ' attention grabber end' ] );
        trial_start_event_id = event_buffer.event_id(trial_start_ebidx);
        trial_gaze_captured_ebtime = event_buffer.time(trial_start_ebidx);
        trial_start_ebtime = event_buffer.time(trial_start_ebidx) - 0.500; 

        
        % audio onset
        trial_audio_ebidx = find(event_buffer.event == ['exp2 trial ', num2str(ttrial), ' sound onset' ] );
        trial_audio_event_id = event_buffer.event_id(trial_audio_ebidx);
        trial_audio_ebtime = event_buffer.time(trial_audio_ebidx);
        
        if isempty(trial_audio_ebidx)
            trial_data.isi(t,1)=NaN;
            trial_data.DurationTrial(t,1)=NaN;
            trial_data.SamplingRate(t,1)=NaN;
            trial_data.DistanceScreen(t,1) = NaN;
            trial_data.MissingRawGazeTrial(t,1) = NaN;
            trial_data.inRSbeforeSound(t,1) = NaN;
            trial_data.ValidInRSbeforeSound(t,1) = 0;
            trial_data.onScreenBeforeSound(t,1) = NaN;
            trial_data.ValidOnScreenBeforeSound(t,1) = 0;
            trial_data.onScreenAfterSound(t,1) = NaN;
            trial_data.ValidOnScreenAfterSound(t,1) = 0;
            trial_data.onScreenAfterSound_p1(t,1) = NaN;
            trial_data.onScreenAfterSound_p2(t,1) = NaN;
            trial_data.ValidOnScreenAfterSound_distributed(t,1) = 0;
            trial_data.DurationArray(t,1) = NaN;
            trial_data.onScreenArray(t,1) = NaN;
            trial_data.onRS(t,1) = NaN;
            trial_data.ValidGaze(t,1) = 0;
            trial_data.DistanceCenter(t,1) = NaN;
            trial_data.PupilEye(t,1) = NaN;
            trial_data.MissingRawPupil(t,1) = NaN;
            trial_data.PupilOutside3STD(t,1) = NaN;
            trial_data.PupilOutsideRange(t,1) = NaN;
            trial_data.MeanBaseline(t,1) = NaN;
            trial_data.MissingBaseline(t,1) = NaN;
            trial_data.MeanResponse(t,1) = NaN;
            trial_data.MissingResponse(t,1) = NaN;
            trial_data.PupilDilation(t,1) = NaN;
            continue
        end

        % visual onset
        trial_visual_ebidx = find(event_buffer.event == ['exp2 trial ', num2str(ttrial), ' video onset' ] );
        trial_visual_event_id = event_buffer.event_id(trial_visual_ebidx);
        trial_visual_ebtime = event_buffer.time(trial_visual_ebidx);
        
        % visual offset = the end of the trial
        trial_end_ebidx = find(event_buffer.event == ['exp2 trial ', num2str(ttrial), ' video offset' ] );
        trial_end_event_id = event_buffer.event_id(trial_end_ebidx);
        trial_end_ebtime = event_buffer.time(trial_end_ebidx);
        
        if isempty(trial_end_ebidx)
            trial_data.isi(t,1)=NaN;
            trial_data.DurationTrial(t,1)=NaN;
            trial_data.SamplingRate(t,1)=NaN;
            trial_data.DistanceScreen(t,1) = NaN;
            trial_data.MissingRawGazeTrial(t,1) = NaN;
            trial_data.inRSbeforeSound(t,1) = NaN;
            trial_data.ValidInRSbeforeSound(t,1) = 0;
            trial_data.onScreenBeforeSound(t,1) = NaN;
            trial_data.ValidOnScreenBeforeSound(t,1) = 0;
            trial_data.onScreenAfterSound(t,1) = NaN;
            trial_data.ValidOnScreenAfterSound(t,1) = 0;
            trial_data.onScreenAfterSound_p1(t,1) = NaN;
            trial_data.onScreenAfterSound_p2(t,1) = NaN;
            trial_data.ValidOnScreenAfterSound_distributed(t,1) = 0;
            trial_data.DurationArray(t,1) = NaN;
            trial_data.onScreenArray(t,1) = NaN;
            trial_data.onRS(t,1) = NaN;
            trial_data.ValidGaze(t,1) = 0;
            trial_data.DistanceCenter(t,1) = NaN;
            trial_data.PupilEye(t,1) = NaN;
            trial_data.MissingRawPupil(t,1) = NaN;
            trial_data.PupilOutside3STD(t,1) = NaN;
            trial_data.PupilOutsideRange(t,1) = NaN;
            trial_data.MeanBaseline(t,1) = NaN;
            trial_data.MissingBaseline(t,1) = NaN;
            trial_data.MeanResponse(t,1) = NaN;
            trial_data.MissingResponse(t,1) = NaN;
            trial_data.PupilDilation(t,1) = NaN;
            continue
        end

        % save the time from audio onset to visual onset which is in
        % interstimuli interval (random between 80-400 ms)
        trial_data.isi(t,1) = trial_audio_ebtime - trial_visual_ebtime;
      
        % Find start and end of trial in indexes on Buffer
        % trial based on 500 ms before gaze captured and visual offset
        %trial_start_gbidx = find(gaze_buffer.event_id >= trial_start_event_id, 1, "first" )
        trial_start_gbidx = find(gaze_buffer.time >= trial_start_ebtime, 1, "first" );
        trial_end_gbidx = find(gaze_buffer.time <= trial_end_ebtime, 1, "last" );
        
        % save total duration of trial - from 500 ms before gaze captured 
        % to visual offset
        trial_data.DurationTrial(t,1) = trial_end_ebtime - trial_start_ebtime;
        
        % segment time buffer
        tb = gaze_buffer.time( trial_start_gbidx:trial_end_gbidx, 1 );
        time_trial = (tb(:) - tb(1))*1000; % convert time buffer to ms from visual onset
        
        % save sampling rate in Hz
        sr = etDetermineSampleRate(tb*1000000); % see function below 
        trial_data.SamplingRate(t,1) = sr;
        
        % segment X buffer
        gaze_lx = gaze_buffer.left_gaze_x( trial_start_gbidx:trial_end_gbidx, 1);
        gaze_rx = gaze_buffer.right_gaze_x( trial_start_gbidx:trial_end_gbidx, 1);
        
        % segment Y buffer
        gaze_ly = gaze_buffer.left_gaze_y( trial_start_gbidx:trial_end_gbidx, 1);
        gaze_ry = gaze_buffer.right_gaze_y( trial_start_gbidx:trial_end_gbidx, 1);  
        
        
        %% Save mean distance to screen
        % Lowe recommended using gaze_buffer.left_eye_cam_y and right to get distance
        distance_ly = gaze_buffer.left_eye_cam_y( trial_start_gbidx:trial_end_gbidx, 1);
        distance_ry = gaze_buffer.right_eye_cam_y( trial_start_gbidx:trial_end_gbidx, 1);
        trial_data.DistanceScreen(t,1) = nanmean( mean( [distance_ly, distance_ry] ,2, 'omitnan') );
      
        
        %% Save missing data before interpolation for the entire trial
        RawMissing = (isnan(gaze_lx) & isnan(gaze_rx)) | (isnan( gaze_ly) & isnan( gaze_ry) );
        trial_data.MissingRawGazeTrial(t,1) = sum(RawMissing) / length(RawMissing) ;
        
        %% blink check

        tmp_gaze_lx = findcontig(isnan(gaze_lx));
        for i=1:length(tmp_gaze_lx(:,1))
            if tmp_gaze_lx(i,3)>75 && tmp_gaze_lx(i,3)<250
                if tmp_gaze_lx(i,1)>25 && tmp_gaze_lx(i,2)<(length(gaze_lx)-25)
                    for ii=tmp_gaze_lx(i,1)-25:tmp_gaze_lx(i,2)+25
                        gaze_lx(ii)=NaN;
                    end
                end
            end
        end

        tmp_gaze_ly = findcontig(isnan(gaze_ly));
        for i=1:length(tmp_gaze_ly(:,1))
            if tmp_gaze_ly(i,3)>75 && tmp_gaze_ly(i,3)<250
                if tmp_gaze_ly(i,1)>25 && tmp_gaze_ly(i,2)<(length(gaze_ly)-25)
                    for ii=tmp_gaze_ly(i,1)-25:tmp_gaze_ly(i,2)+25
                        gaze_ly(ii)=NaN;
                    end
                end
            end
        end

        tmp_gaze_rx = findcontig(isnan(gaze_rx));
        for i=1:length(tmp_gaze_rx(:,1))
            if tmp_gaze_rx(i,3)>75 && tmp_gaze_rx(i,3)<250
                if tmp_gaze_rx(i,1)>25 && tmp_gaze_rx(i,2)<(length(gaze_rx)-25)
                    for ii=tmp_gaze_rx(i,1)-25:tmp_gaze_rx(i,2)+25
                        gaze_rx(ii)=NaN;
                    end
                end
            end
        end

        tmp_gaze_ry = findcontig(isnan(gaze_ry));
        for i=1:length(tmp_gaze_ry(:,1))
            if tmp_gaze_ry(i,3)>75 && tmp_gaze_ry(i,3)<250
                if tmp_gaze_ry(i,1)>25 && tmp_gaze_ry(i,2)<(length(gaze_ry)-25)
                    for ii=tmp_gaze_ry(i,1)-25:tmp_gaze_ry(i,2)+25
                        gaze_ry(ii)=NaN;
                    end
                end
            end
        end


        %% interpolate and average eyes
        % Data was interpolated linearly over gaps in the data shorter 
        % than 150 ms (as in Kleberg 2019, same as CBCD/Luke Mason's scripts)
        
        [mb, flags] = etInterpBuffer(gaze_lx,gaze_ly, gaze_rx,gaze_ry, tb, interp_maxGap); % see function below
        
        lx_out = mb(:, 1);
        ly_out = mb(:, 2);
        rx_out = mb(:, 3);
        ry_out = mb(:, 4);
        
        gaze_x_trial = mean( [lx_out, rx_out] ,2, 'omitnan');
        gaze_y_trial = mean( [ly_out, ry_out] ,2, 'omitnan');
        
%         % get missing data after interpolation
%         gaze_missing = isnan(gaze_x_trial) | isnan( gaze_y_trial) ;
        

        %% Validate trial
        
        % compute binary vectors for gaze on CS, on screen, outside screen
        
        inAOIfixation =...
            gaze_x_trial >= fixation_AOI(1) &...
            gaze_x_trial <= fixation_AOI(3) &...
            gaze_y_trial >= fixation_AOI(2) &...
            gaze_y_trial <= fixation_AOI(4);
        
        onScreen = ...
            gaze_x_trial <= width_half &...
            gaze_x_trial >= -width_half &...
            gaze_y_trial <= height_half &...
            gaze_y_trial >= -height_half;
       
        onRightScreen = ...
            gaze_x_trial <= width_half &...
            gaze_x_trial > 0 &...
            gaze_y_trial <= height_half &...
            gaze_y_trial >= -height_half;
        
        outScreen = ...
            gaze_x_trial > width_half |...
            gaze_x_trial < -width_half |...
            gaze_y_trial > height_half |...
            gaze_y_trial < -height_half;
        
        % 1. Trial invalid if interpolated gaze was not at the right side
        % of the screen for at least 40% of the 100 ms before and 500ms
        % after sound onset        
        sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
        check_rs_before = 100; %ms
        check_rs_after = 500;
        crit_minInRS = .4;
        InRS_idx = sound_start_gbidx - round(check_rs_before*sr/1000) : sound_start_gbidx + round(check_rs_after*sr/1000);
        trial_data.inRSbeforeSound(t,1) = sum( onRightScreen(InRS_idx))/ length(onRightScreen(InRS_idx)) ;
        trial_data.ValidInRSbeforeSound(t,1) = trial_data.inRSbeforeSound(t,1) >= crit_minInRS;
         
              
        % 2. Trial invalid if interpolated gaze was missing or outside the screen
        % (ie it was not inside the screen) during the 500 ms before sound
        % onset and for the rest of the video (3000ms)     
        check_OnScreen_before = 500; %ms
        check_OnScreen_after = 0;
        crit_minOnScreen = 0.7;
        OnScreen_idx = sound_start_gbidx - round(check_OnScreen_before*sr/1000) : sound_start_gbidx+ round(check_OnScreen_after*sr/1000); 
        trial_data.onScreenBeforeSound(t,1) = sum( onScreen(OnScreen_idx)) / length( onScreen(OnScreen_idx));
        trial_data.ValidOnScreenBeforeSound(t,1) = trial_data.onScreenBeforeSound(t,1) >= crit_minOnScreen;

        % 3. Trial invalid if interpolated gaze was missing or outside the screen
        % (ie it was not inside the screen) during the 3000ms after sound onset   
        % NB now adapted to late stage time-window (second half
        % post-sound), for the main analyses, OnScreen_before = 0 &
        % OnScreen_after = 2890
        check_OnScreen_before = 1500; %ms before 0 in the main analysis!!!
        check_OnScreen_after = 2890; % 3000;
        crit_minOnScreen = 0.5;
        OnScreen_idx = sound_start_gbidx + round(check_OnScreen_before*sr/1000) : sound_start_gbidx+ round(check_OnScreen_after*sr/1000); %length(onScreen); 
        trial_data.onScreenAfterSound(t,1) = sum( onScreen(OnScreen_idx)) / length( onScreen(OnScreen_idx));
        trial_data.ValidOnScreenAfterSound(t,1) = trial_data.onScreenAfterSound(t,1) >= crit_minOnScreen;

        % 4. Trial invalid if interpolated gaze was missing or outside the screen
        % (ie it was not inside the screen) disproportionately during the 3000ms after sound
          
        check_OnScreen_before = 0; %ms
        check_OnScreen_middle = 1490;
        check_OnScreen_after = 2890; % 3000;
        crit_minOnScreen_p1 = 0.4;
        OnScreen_idx_pre = sound_start_gbidx - round(check_OnScreen_before*sr/1000) : sound_start_gbidx+ round(check_OnScreen_middle*sr/1000); %length(onScreen);
        OnScreen_idx_post = sound_start_gbidx + round(check_OnScreen_middle*sr/1000) : sound_start_gbidx+ round(check_OnScreen_after*sr/1000); %length(onScreen); 
        
        trial_data.onScreenAfterSound_p1(t,1) = sum( onScreen(OnScreen_idx_pre)) / length( onScreen(OnScreen_idx));
        trial_data.onScreenAfterSound_p2(t,1) = sum( onScreen(OnScreen_idx_post)) / length( onScreen(OnScreen_idx));
        
        trial_data.ValidOnScreenAfterSound_distributed(t,1) = trial_data.onScreenAfterSound_p2(t,1)-0.15 < trial_data.onScreenAfterSound_p1(t,1) || trial_data.onScreenAfterSound_p1(t,1) < trial_data.onScreenAfterSound_p2(t,1)+0.15;
        
                    
        
        %% segment visual array presentation
        % from visual onset to visual offset
        
        visual_start_gbidx = find(tb >= trial_visual_ebtime, 1, "first" );
        visual_end_gbidx = find(tb <= trial_end_ebtime, 1, "last" );
        
        idx_visual = visual_start_gbidx:visual_end_gbidx;
        
        gaze_x = gaze_x_trial(idx_visual);
        gaze_y = gaze_y_trial(idx_visual);
        time = tb(idx_visual);
        
        % save total duration of presentation in seconds from number of
        % samples collected and sampling rate (= 4 seconds)
        trial_data.DurationArray(t,1) = size(gaze_x, 1) / sr;
        
%         % get missing data after interpolation during the presentation
%         missing = isnan(gaze_x) | isnan( gaze_y);
                

        %% VISUAL RESPONSE STATS
        
               
        onScreen = ...
            gaze_x <= width_half &...
            gaze_x >= -width_half &...
            gaze_y <= height_half &...
            gaze_y >= -height_half;
        
        onRightScreen = ...
            gaze_x <= width_half &...
            gaze_x > 0 &...
            gaze_y <= height_half &...
            gaze_y >= -height_half;

        % save total time on Screen from number of samples collected and
        % sampling rate (in seconds)
        trial_data.onScreenArray(t,1) = sum( onScreen) / sr;
             
        % get time on each AOI from number of samples collected and
        % sampling rate (in seconds)
        trial_data.onRS(t,1) = sum(onRightScreen) / sr;
        
                       
        %% Validate trial based on all flags
        trial_data.ValidGaze(t,1) = trial_data.ValidInRSbeforeSound(t,1) == 1 &&...
            trial_data.ValidOnScreenBeforeSound(t,1) == 1 &&...
            trial_data.ValidOnScreenAfterSound(t,1) == 1 &&...
            trial_data.ValidOnScreenAfterSound_distributed(t,1) == 1;       
        
        trial_data.DistanceCenter(t,1) = mean(sqrt(gaze_x.^2 + gaze_y.^2),1,'omitnan');
        %% PUPIL RESPONSE STATS
        
        % segment pupil data
        pupilL = gaze_buffer.left_pupil_measure1( trial_start_gbidx:trial_end_gbidx, 1);
        pupilR = gaze_buffer.right_pupil_measure1( trial_start_gbidx:trial_end_gbidx, 1);
        
        % missing pupil
        pupilL_miss = isnan(pupilL);
        pupilR_miss = isnan(pupilR);
        
        % smooth pupil data with a moving window average (100 ms window)
        % using matlab function movmean(, 'omitnan') but replacing all 
        % missing data with NaN again after        
        mawindow_ms = 100; % 100 ms 
        window = (sr*mawindow_ms) / 1000; % in samples
        pupilL_smooth = movmean(pupilL,window, 'omitnan');
        pupilL_smooth(pupilL_miss) = NaN;
        pupilR_smooth = movmean(pupilR,window, 'omitnan');
        pupilR_smooth(pupilR_miss) = NaN;
        
        % average the pupil signal 
        % if we have data for both eyes
        if sum(~pupilL_miss) > 0 &&  sum(~pupilR_miss) > 0 && sum(~pupilL_miss & ~pupilR_miss) > 0
            
            % average the pupil signal considering a dynamic offset
            % using a function from pupil-size-master (Kret, & Sjak-Shie, 2019)
            pupil_mean = genMeanDiaSamples(time_trial, pupilL_smooth, pupilR_smooth, ~pupilL_miss, ~pupilR_miss);

            % Ana found some errors in this function when the begining of 
            % the signal was missing / only one eye was present so in some
            % cases the normal mean was computed instead
            
            % check which data stream has more data
            streams = ['L', 'R', 'Mean'];
            streams_miss = [sum(pupilL_miss), sum(pupilR_miss), sum(isnan(pupil_mean)) ];
            idx_stream = find( streams_miss == min(streams_miss) );
            
            if any(idx_stream == 3) && ~isempty(pupil_mean)
                trial_data.PupilEye(t,1) = 1; % we used mean pupil with dynamic offset
            else
                pupil_mean =  mean( [pupilL_smooth, pupilR_smooth] ,2, 'omitnan');
                trial_data.PupilEye(t,1) = 2; % we used normal mean
            end
            
        else
            % if one eye is completely missing do normal average which is going
            % to take only one eye.
            pupil_mean =  mean( [pupilL_smooth, pupilR_smooth] ,2, 'omitnan');
            trial_data.PupilEye(t,1) = 3; % we used only one eye
            
            %         % check which data stream as more data
            %         streams = ['L', 'R'];
            %         streams_miss = [sum(pupilL_miss) ,  sum(pupilR_miss)];
            %
            %         idx_stream = find( streams_miss == min(streams_miss) );
            %
            %         if idx_stream == 1
            %             trial_data.Pupil_eye(t,1) = 'Only L     '; % we used left pupil
            %         elseif idx_stream == 2
            %             trial_data.Pupil_eye(t,1) = 'Only R     '; % we used right pupil
            %         end
        end
        
        % save how much missing data before interpolation we had
        trial_data.MissingRawPupil(t,1) = sum(isnan(pupil_mean)) / length(pupil_mean) ;
        
        % Exclude invalid samples based on reasonable mm range (Kret, & Sjak-Shie, 2019)
        pupil_mean(pupil_mean < 1.5) = NaN;
        pupil_mean(pupil_mean > 9) = NaN;
        
        % Exclude invalid samples based on outliers (Mathôt & Vilotijević, 2022)
        Pmean = double( nanmean(pupil_mean));
        P3STD = double( 3*std(pupil_mean, 'omitnan'));
        pupil_mean(pupil_mean > Pmean + P3STD) = NaN;
        pupil_mean(pupil_mean < Pmean - P3STD) = NaN;
        
        % Save how much invalid data we excluded
        trial_data.PupilOutside3STD(t,1) = sum( pupil_mean > Pmean + P3STD | pupil_mean < Pmean - P3STD ) / length(pupil_mean);
        trial_data.PupilOutsideRange(t,1) = sum( pupil_mean < 1.5 | pupil_mean > 9 ) / length(pupil_mean);
        
        % Exclude outside screen samples
        pupil_mean(outScreen) = NaN;
        
        
        % Interpolate missing/invalid pupil data
        % linearly over gaps in the data shorter 
        % than 150 ms (Kleberg 2019, same as CBCD/Luke Mason's scripts)
        pupil_miss = isnan(pupil_mean);
        
        gaps = findcontig(pupil_miss);%, true);
        
        % check there are missing gaps and also some pupil data
        if size(gaps,1) > 0 && size(time_trial(~pupil_miss), 1) > 2 && ~isempty(interp_maxGap)
            
            pupil_interp = interp1(time_trial(~pupil_miss) ...
                ,pupil_mean(~pupil_miss)...
                ,time_trial,'linear');
            
            % select those runs that are higher than the maximum length that we
            % want to interpolate over
            large_gaps = gaps( gaps(:,3)* (1000 / sr) > interp_maxGap, : );
            
            % loop through all invalid interpolated sections (i.e. missing
            % data was more than criterion) and get where they happened
            indexes_large_gaps = [];
            for i = 1:size(large_gaps,1)
                indexes_large_gaps = horzcat(indexes_large_gaps, linspace(large_gaps(i,1),large_gaps(i,2),large_gaps(i,3)));
            end 
            pupil_interp(indexes_large_gaps) = NaN;  % Set all values in large gaps to NaN again
        else
            
            % if there are no missing gaps just keep the mean pupil
            pupil_interp = pupil_mean;
        end
        
        
        % BASELINE AND RESPONSE STATS
        
        % save mean pupil size in the baseline
        % as well as missing after interpolation
        % baseline is defined as the 500 ms before audio onset
        sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
        duration_baseline_before = 500; %ms
        duration_baseline_after = 0;
        baseline_idx = sound_start_gbidx - round(duration_baseline_before*sr/1000) : sound_start_gbidx + round(duration_baseline_after*sr/1000);
        trial_data.MeanBaseline(t,1) = mean( pupil_interp(baseline_idx)  , 'omitnan');
        trial_data.MissingBaseline(t,1) = sum( isnan( pupil_interp(baseline_idx)) ) / size(baseline_idx,2);
      
% %         % save mean pupil size in the response
% %         % as well as missing after interpolation
% %         % response is defined as the 3 seconds after sound onset
% %         sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
% %         duration_response_before = 0; %ms
% %         duration_response_after = 1500; % WAS 3000MS BEFORE final 2980!!!
% %         response_idx = sound_start_gbidx - round(duration_response_before*sr/1000) : sound_start_gbidx + round(duration_response_after*sr/1000); %length(pupil_interp);
% %         trial_data.MeanResponse(t,1) = mean( pupil_interp(response_idx)  , 'omitnan');
% %         trial_data.MissingResponse(t,1) = sum( isnan( pupil_interp(response_idx)) ) / size(response_idx,2);

        % SECOND TIME-WINDOW save mean pupil size in the response
        % as well as missing after interpolation
        % response is defined as the 3 seconds after sound onset
        sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
        duration_response_before = 1500; %ms
        duration_response_after = 2890; % WAS 3000MS BEFORE final 2980!!!
        response_idx = sound_start_gbidx + round(duration_response_before*sr/1000) : sound_start_gbidx + round(duration_response_after*sr/1000); %length(pupil_interp);
        trial_data.MeanResponse(t,1) = mean( pupil_interp(response_idx)  , 'omitnan');
        trial_data.MissingResponse(t,1) = sum( isnan( pupil_interp(response_idx)) ) / size(response_idx,2);
        
        % save pupil dilation indexed as mean pupil response - mean pupil baseline
        if ~isnan(trial_data.MeanBaseline(t,1)) && ~isnan(trial_data.MeanResponse(t,1))
            trial_data.PupilDilation(t,1) = trial_data.MeanResponse(t,1) - trial_data.MeanBaseline(t,1);
        else
            % if we are missing baseline or response mark as missing
            trial_data.PupilDilation(t,1) = NaN;
        end
        
        % Collect normalized pupil during baseline and during response for
        % later average tracing / plots
        pupil_normal = (pupil_interp - trial_data.MeanBaseline(t,1)) / trial_data.MeanBaseline(t,1);
        
        if t_for_B==1
           pupil_all.response(:, t_for_B) = pupil_normal(response_idx);
           pupil_all.baseline(:, t_for_R) = pupil_normal(baseline_idx);
        end
        if t_for_B>1 && size(pupil_all.response(:, t_for_B-1),1)==size(pupil_normal(response_idx),1)
           pupil_all.response(:, t_for_B) = pupil_normal(response_idx);
           pupil_all.baseline(:, t_for_R) = pupil_normal(baseline_idx);
        elseif t_for_B>1 && size(pupil_all.response(:, t_for_B-1),1)~=size(pupil_normal(response_idx),1)
           pupil_all.response(:, t_for_B) = NaN(size(pupil_all.response(:, t_for_B-1),1),1);
           pupil_all.baseline(:, t_for_R) = NaN(size(pupil_all.baseline(:, t_for_R-1),1),1);
        end
        t_for_B = t_for_B +1;
        t_for_R = t_for_R +1;
 
        %% plot trial data
        if plot_figure &&...
                any(t == trials_to_plot)

            % plot gaze data
            spx = subplot(2, 2, 1);
            hold(spx, 'on')
            title('Gaze plot')
            plot(tb, gaze_lx, '-r', 'linewidth', 2)
            plot(tb, gaze_ly, '-g', 'linewidth', 2)
            
            plot(tb, gaze_rx, '-r', 'linewidth', 2)
            plot(tb, gaze_ry, '-g', 'linewidth', 2)
            
            plot(tb, gaze_x_trial, '-k', 'linewidth', 1)
            plot(tb, gaze_y_trial, '-k', 'linewidth', 1)
            
            % plot inCS validation period before sound onset
            sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
            check_rs_before = 500; %ms
            check_rs_after = 0;
            pos_sound = [trial_audio_ebtime - (check_rs_before/1000),...
                -width_half,...
                (check_rs_before/1000 + check_rs_after/1000),...
                width_half*2];
            rectangle('Position',pos_sound, 'FaceColor',[0.4940 0.1840 0.5560 0.20])
            
            % plot inCS/ onScreen validation period after onset
            check_rs_before = 0; %ms
            check_rs_after = 3000;
            pos_visual = [trial_audio_ebtime - (check_rs_before/1000),...
                -width_half,...
                (check_rs_before/1000 + check_rs_after/1000),...
                width_half*2];
            rectangle('Position',pos_visual, 'FaceColor',[0.3010 0.7450 0.9330 0.20])
            
            % plot sound onset
            line([trial_audio_ebtime, trial_audio_ebtime],[-width_half, width_half], 'color', [0, 0, 0],...
                'linewidth', 2);
            % plot visual onset
            line([trial_visual_ebtime, trial_visual_ebtime],[-width_half, width_half], 'color', [0, 0, 0],...
                'linewidth', 2);
            
            legend({'X gaze','Y gaze'},'Location','southwest')
            
            ylim([-width_half, width_half]) % coordinates are X screen ones
            xlim([tb(1), tb(end)])
            
            % another plot: plot inAOI vectors
            spx = subplot(2, 2, 3);
            hold(spx, 'on')
            title('AOI Hit plot')
            
            scatter(time, onRightScreen*2, 100, 'g', 'filled') 
            scatter(tb, inAOIfixation*1, 100, 'r', 'filled')
            
            line([trial_audio_ebtime, trial_audio_ebtime], [0.5,5.5], 'color', [0, 0, 0],...
                'linewidth', 2);
            line([trial_visual_ebtime, trial_visual_ebtime], [0.5,5.5], 'color', [0, 0, 0],...
                'linewidth', 2);
            
            ylim([0.5,5.5])
            xlim([tb(1), tb(end)])
            
            legend({'Right screen', 'Fixation'},'Location','northwest')

            % another plot: plot pupil tracing
            spx = subplot(2, 2, 2);
            hold(spx, 'on')
            title('Pupil size plot')
            
            plot(tb, pupilL, '-b', 'linewidth', 2)
            plot(tb, pupilR, '-m', 'linewidth', 2)
            plot(tb, pupil_interp, '-k', 'linewidth', 1)
            
            min_pupil = min([pupilL;pupilR; pupil_interp]);
            max_pupil = max([pupilL;pupilR; pupil_interp]);
            
            if ~isnan(min_pupil)
                % pos = [x y w h]
                % plot baseline period
                pos_baseline = [trial_audio_ebtime - (duration_baseline_before/1000),...
                    min_pupil,...
                    (duration_baseline_before/1000 + duration_baseline_after/1000),...
                    max_pupil - min_pupil];
                rectangle('Position',pos_baseline, 'FaceColor',[0.9290 0.6940 0.1250 0.20])
                
                % plot response period
                pos_response = [trial_audio_ebtime - (duration_response_before/1000),...
                    min_pupil,...
                    (duration_response_before/1000 + duration_response_after/1000),...
                    max_pupil - min_pupil];
                rectangle('Position',pos_response, 'FaceColor',[0.8500 0.3250 0.0980 0.20])
                
                line([trial_audio_ebtime, trial_audio_ebtime], [min_pupil-.1, max_pupil+.1], 'color', [0, 0, 0],...
                    'linewidth', 2);
                line([trial_visual_ebtime, trial_visual_ebtime], [min_pupil-.1, max_pupil+.1], 'color', [0, 0, 0],...
                    'linewidth', 2);
                xlim([tb(1), tb(end)])
                ylim([min_pupil-.1, max_pupil+.1])
            end
            
            % another plot: plot some stats
            spx = subplot(2, 2, 4);
            % write down stats
            hold(spx, 'on')
            title('Stats')

            text(0.5, 0.9, ['SR =', num2str( trial_data.SamplingRate(t,1) ), 'Hz'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.8, ['% In right side of screen before Sound onset =', num2str( round( trial_data.inRSbeforeSound(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.8, ['% On Screen after Sound onset =', num2str( round( trial_data.onScreenAfterSound(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.1, 0.7, ['% On Screen after Sound onset p1 =', num2str( round( trial_data.onScreenAfterSound_p1(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.7, ['% On Screen after Sound onset p2 =', num2str( round( trial_data.onScreenAfterSound_p2(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.5, ['% Miss Pupil during Baseline =', num2str( round( trial_data.MissingBaseline(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.5, ['% Miss Pupil during Response =', num2str( round( trial_data.MissingResponse(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            
            % mean and STD of pupil size, max and min
            text(0.1, 0.4, ['Pupil Mean =', num2str(mean(pupil_interp, 'omitnan'))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.4, ['Pupil STD =', num2str(std(pupil_interp, 'omitnan'))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.1, 0.3, ['Pupil Minimum =', num2str(min(pupil_interp))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.3, ['Pupil Maximum =', num2str(max(pupil_interp))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.2, ['Pupil outside 3STD =', num2str( round( trial_data.PupilOutside3STD(t,1) *100 )), '%'],...
                'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.2, ['Pupil outside Range (1.5-9) =', num2str( round( trial_data.PupilOutsideRange(t,1)  *100 )), '%'],...
                'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.3, 0.1, ['Valid Trial? ', num2str(trial_data.ValidGaze(t,1)), ' (1=Yes, 0=No)'],...
                'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            
            set(gcf, 'WindowState', 'fullscreen')
            saveas(gcf,['../results_habituation_shortwindow2/', id_name, ' trial ', num2str(ttrial),' plots.png'])
            %   pause;
            close(gcf)
            
            
            % another figure: plot "heatmap"
            figure, scatter(gaze_x, gaze_y)
            hold on
            rectangle('Position',[-stim_box_grab/2 -stim_box_grab/2 stim_box_grab/2 stim_box_grab/2]) % fixation AOI
            rectangle('Position',[0, -height_half, width_half, 2*height_half]) % habituation AOI
            
            
            xlim([-width_half, width_half])
            ylim([-height_half, height_half])
            
            set(gcf, 'WindowState', 'fullscreen')
            saveas(gcf,['../results_habituation_shortwindow2/', id_name, ' trial ', num2str(ttrial),' heatmap.png'])
            %pause;
            close(gcf)
            close all            
        end        
    end
    
    % % uncomment to save trial data in an invidiual file
    % if save_data == 1
    %   writetable(struct2table(trial_data), strcat(id_name, '_', experiment_date, '_Dataset_long.csv') )
    % end
    
    %% collect all participants data in a struct

    if ~exist('trial_data_all') % on first subject
        trial_data_all = trial_data;
    else
        % merge trial data for subject with the rest of the trial data
        trial_fields = fieldnames(trial_data);
        for i = 1:length(trial_fields)
            trial_data_all.(trial_fields{i}) = [trial_data_all.(trial_fields{i}); trial_data.(trial_fields{i})];
        end
    end
    
    %% collect average stats for each individual (only valid trials!)
        % see variables saved described below
        
        % select valid trials
        idx_Valid = trial_data.ValidGaze == 1;
        
        % save total number of trials
        trial_data_all_wide.NumberTrials(f,1) = length(trial_data.id_name);
        
        trial_data_all_wide.MissingRawGazeTrial(f,1) = mean(trial_data.MissingRawGazeTrial(idx_Valid));
        trial_data_all_wide.onScreenArray(f,1) = mean(trial_data.onScreenArray(idx_Valid));
        
        outcomes = {'PupilDilation'};
        
        % select trials for each condition 
        clear idx
        idx.multi1 = trial_data.StimulusType == 1  & trial_data.Order == 1 ;
        idx.multi2 = trial_data.StimulusType == 1  & trial_data.Order == 2 ;
        idx.multi3 = trial_data.StimulusType == 1  & trial_data.Order == 3 ;
        idx.multi4 = trial_data.StimulusType == 1  & trial_data.Order == 4 ;
        idx.multi5 = trial_data.StimulusType == 1  & trial_data.Order == 5 ;
        idx.multi6 = trial_data.StimulusType == 1  & trial_data.Order == 6 ;

        idx.uni1 = trial_data.StimulusType == 2  & trial_data.Order == 1 ;
        idx.uni2 = trial_data.StimulusType == 2  & trial_data.Order == 2 ;
        idx.uni3 = trial_data.StimulusType == 2  & trial_data.Order == 3 ;
        idx.uni4 = trial_data.StimulusType == 2  & trial_data.Order == 4 ;
        idx.uni5 = trial_data.StimulusType == 2  & trial_data.Order == 5 ;
        idx.uni6 = trial_data.StimulusType == 2  & trial_data.Order == 6 ;
        
        
        f_idx = fieldnames(idx); % each possible combination of conditions
    
        % loop through all possible combination of conditions and save number
        % of valid trials
        for c=1:size(f_idx,1) % 10 Ns 
            trial_data_all_wide.(['NumberValidTrials', f_idx{c}])(f,1) = sum(idx.(f_idx{c}) & idx_Valid);
        end
        
        % loop through all measures and possible combination of conditions and 
        % save mean across valid trials
        for o = 1:size(outcomes,1) % 1 outcome here: pupil dilation
            for c=1:size(f_idx,1) % 12 conditions = 12 means
                if trial_data_all_wide.(['NumberValidTrials', f_idx{c}])(f,1) > 0
                    trial_data_all_wide.(['Mean',outcomes{o}, f_idx{c}])(f,1) = nanmean( trial_data.(outcomes{o})(idx.(f_idx{c}) &...
                        idx_Valid) );
                else
                    % if there are no valid trials mark as missing
                    trial_data_all_wide.(['Mean',outcomes{o}, f_idx{c}])(f,1) = NaN;
                end
            end
        end
    
        % loop through all measures except FirstLookFace and possible 
        % combination of conditions and save std across valid trials
        for o = 1:size(outcomes,1) % 1 outcome here: pupil dilation
            for c=1:size(f_idx,1) % 12 conditions = 12 stds
                if trial_data_all_wide.(['NumberValidTrials', f_idx{c}])(f,1) > 0
                    trial_data_all_wide.(['STD',outcomes{o}, f_idx{c}])(f,1) = nanstd( trial_data.(outcomes{o})(idx.(f_idx{c}) &...
                        idx_Valid) );
                else
                    % if there are no valid trials mark as missing
                    trial_data_all_wide.(['STD',outcomes{o}, f_idx{c}])(f,1) = NaN;
                end
            end
        end

    
    % clear vars
    clear trial_data
    clearvars -except trial_data_all... 
        id interp_maxGap multi uni width_half height_half plot_figure save_data stim_box_grab stim_box...
        resultsPath resultsPath_calib files...
        pupil_all t_for_B t_for_R
%trial_data_all_wide...
    
end

if save_data == 1
    
    %% write long file
    writetable(struct2table(trial_data_all), ['../results_habituation_shortwindow2/', date, ' All_Datasets_long.csv'] )
    
    % % header are:
    % id_name	Participant ID
    % id_date	Experiment date
    % id_trial	Trial nr
    % StimulusType	1 = multi-sensory event, 2 = uni-sensory event
    % DurationTrial	Time (seconds) from 500 ms before sound onset was triggered (gaze was captured) to visual array offset
    % SamplingRate	Number of samples per second, 600 Hz
    % MissingRawGazeTrial	Proportion raw (before interpolation) missing gaze data in the trial (from 500 ms before sound onset to visual array offset)
    % inRSbeforeSound	Trial validity flag 1: proportion of interpolated gaze that was at the right side of the screen during the 500 ms before audio onset
    % ValidInRSbeforeSound	Trial validity flag 1: 0 = invalid = inRSbeforeSound < .4, 1 = valid
    % onScreenBeforeSound	Trial validity flag 3: proportion of interpolated gaze that was inside the screen during the 500 ms before sound onset
    % ValidOnScreenBeforeSound	Trial validity flag 3: 0 = invalid = onScreenAfterSound < 1, 1 = valid
    % onScreenAfterSound	Trial validity flag 3: proportion of interpolated gaze that was inside the screen during the 3000ms after sound onset
    % ValidOnScreenAfterSound	Trial validity flag 3: 0 = invalid = onScreenAfterSound < 1, 1 = valid
    % onScreenAfterSound_p1	Trial validity flag 4: proportion of interpolated gaze that was inside the screen during the first 1500ms after sound onset
    % onScreenAfterSound_p2	Trial validity flag 4: proportion of interpolated gaze that was inside the screen during the second 1500ms after sound onset
    % ValidOnScreenAfterSound_distributed	Trial validity flag 4: 0 =
    % invalid, 1 = valid:
    % onScreenAfterSound_p2-0.15<onScreenAfterSound_p1<onScreenAfterSound_p2+0.15
    % DurationArray	Time (seconds) from visual onset to visual array offset
    % onScreenArray	Cumulative time (seconds) in the screen during visual presentation (visual onset to visual offset)
    % ValidGaze	Whether the trial is valid: 0 = invalid = any of the above flags are 0, 1 = valid
    % PupilEye	1 = mean pupil with dynamic offset, 2 = normal mean pupil, 2 = one eye pupil
    % MissingRawPupil	Proportion raw (before interpolation) missing pupil data in the trial (from 500 ms before sound onset to visual array offset)
    % PupilOutside3STD	Proportion of pupil samples that were excluded because they were outliers (above or below 3STD of trial mean)
    % PupilOutsideRange	Proportion of pupil samples that were excluded because they were invalid (below or above 1.5-9 mm)
    % MeanBaseline	Mean pupil size (mm) during baseline period
    % MissingBaseline	Proportion missing data after interpolation during baseline
    % MeanResponse	Mean pupil size during (mm) response period
    % MissingResponse	Proportion missing data after interpolation during response
    % PupilDilation	Mean Baseline from Mean Response (mm)
    
    %% write wide file
    %writetable(struct2table(trial_data_all_wide), ['../results_habituation/', date, ' All_Datasets_wide.csv'] )
    
    % % headers are:
    % subject	Participant ID
    % id_name	Participant ID with time
    % calib_time	Calibration date and time
    % mean_distance_gaze_to_target	Mean distance (pixels) from gaze to target during calibration validation
    % std_distance_gaze_to_target	STD distance (pixels) from gaze to target during calibration validation
    % proportion_of_time_gaze_on_screen	Proportion gaze on the screen during calibration validation
    % NumberTrials	Number of trials presented
    % MissingRawGazeTrial	Mean proportion raw (before interpolation) missing gaze data across valid trials (from 500 ms before sound onset to visual array offset)
    % onScreenArray	Mean cumulative time (seconds) in the screen during visual presentation (visual onset to visual offset) across valid trials
    % NumberValidTrialsSilent1	Number of valid trials in condition
    % NumberValidTrialsSilent2	Number of valid trials in condition
    % NumberValidTrialsSocialLow1	Number of valid trials in condition
    % NumberValidTrialsSocialLow2	Number of valid trials in condition
    % NumberValidTrialsSocialHigh1	Number of valid trials in condition
    % NumberValidTrialsSocialHigh2	Number of valid trials in condition
    % NumberValidTrialsNonSocialLow1	Number of valid trials in condition
    % NumberValidTrialsNonSocialLow2	Number of valid trials in condition
    % NumberValidTrialsNonSocialHigh1	Number of valid trials in condition
    % NumberValidTrialsNonSocialHigh2	Number of valid trials in condition
    % MeanMinLatencySilent1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySilent2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialLow1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialLow2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialHigh1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialHigh2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialLow1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialLow2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialHigh1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialHigh2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanFirstLookFaceSilent1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSilent2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialLow1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialLow2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialHigh1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialHigh2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialLow1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialLow2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialHigh1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialHigh2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanPupilDilationSilent1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSilent2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialLow1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialLow2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialHigh1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialHigh2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialLow1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialLow2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialHigh1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialHigh2	Mean pupil dilation (response-baseline) across valid trials in condition
    % STDMinLatencySilent1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySilent2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialLow1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialLow2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialHigh1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialHigh2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialLow1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialLow2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialHigh1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialHigh2	STD latency (seconds) of first look detected across valid trials in condition
    % STDPupilDilationSilent1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSilent2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialLow1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialLow2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialHigh1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialHigh2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialLow1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialLow2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialHigh1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialHigh2	STD pupil dilation (response-baseline) across valid trials in condition

end

%% plot average pupil trace

if 0
    
    pupil_baseline = mean( pupil_all.baseline ,2, 'omitnan');
    tb_baseline = linspace(-280,-80,size(pupil_baseline, 1)); % edit if baseline timing changes
    
    pupil_response = mean( pupil_all.response ,2, 'omitnan');
    tb_response = linspace(0,2000,size(pupil_response, 1)); % edit if response timing changes
    
    min_pupil = min([pupil_baseline;pupil_response]);
    max_pupil = max([pupil_baseline;pupil_response]);
    
    % this needs to be adapted depending on the conditions that are
    % relevant
    % Elin/ Matilda requested social vs non-social and repetition effect        
    idx_multi1 = trial_data_all.StimulusType == 1  & trial_data_all.Order == 1 ;
    idx_multi2 = trial_data_all.StimulusType == 1  & trial_data_all.Order == 2 ;
    idx_multi3 = trial_data_all.StimulusType == 1  & trial_data_all.Order == 3 ;
    idx_multi4 = trial_data_all.StimulusType == 1  & trial_data_all.Order == 4 ;
    idx_multi5 = trial_data_all.StimulusType == 1  & trial_data_all.Order == 5 ;
    idx_multi6 = trial_data_all.StimulusType == 1  & trial_data_all.Order == 6 ;

    idx_uni1 = trial_data_all.StimulusType == 2  & trial_data_all.Order == 1 ;
    idx_uni2 = trial_data_all.StimulusType == 2  & trial_data_all.Order == 2 ;
    idx_uni3 = trial_data_all.StimulusType == 2  & trial_data_all.Order == 3 ;
    idx_uni4 = trial_data_all.StimulusType == 2  & trial_data_all.Order == 4 ;
    idx_uni5 = trial_data_all.StimulusType == 2  & trial_data_all.Order == 5 ;
    idx_uni6 = trial_data_all.StimulusType == 2  & trial_data_all.Order == 6 ;
 
    
    idx = [idx_multi1, idx_multi2, idx_multi3, idx_multi4, idx_multi5, idx_multi6, idx_uni1, idx_uni2, idx_uni3, idx_uni4, idx_uni5, idx_uni6];
    idx_color = ['-r'; '-b'; '-g'; '-m'; '-y'; '-k'; ':r'; ':b'; ':g'; ':m'; ':y'; ':k'];

%% loop through conditions defined in idx 

for c = 1:12
    pupil_baseline_mean = mean( pupil_all.baseline(:, idx(:,c)) ,2, 'omitnan');
    pupil_baseline_std =  std( pupil_all.baseline(:, idx(:,c)),0, 2 , 'omitnan') / sqrt(length(pupil_all.baseline(:, idx(:,c))));
    
    % plot the average trace for the specific condition
    plot(tb_baseline, pupil_baseline_mean,idx_color(c,:), 'linewidth', 2)
    hold on 
end

for c = 1:12
    pupil_response_mean = mean( pupil_all.response(:, idx(:,c)) ,2, 'omitnan');
    pupil_response_std =  std( pupil_all.response(:, idx(:,c)),0, 2 , 'omitnan') / sqrt(length(pupil_all.response(:, idx(:,c))));
    
    plot(tb_response, pupil_response_mean, idx_color(c,:), 'linewidth', 2)
    hold on
end

%     line([-80, -80], [-.1, 0.1], 'color', [0, 0, 0],...
%         'linewidth', 2);
    line([0, 0], [-.1, 0.1], 'color', [0, 0, 0],...
        'linewidth', 2);
    
    pos_baseline = [-500,-.08, 0,.1];
    rectangle('Position',pos_baseline, 'FaceColor',[0.8 0.8 0.8 0.25])
    
    legend('Multisensory 1', 'Multisensory 2', 'Multisensory 3', 'Multisensory 4','Multisensory 5','Multisensory 6','Unisensory 1', 'Unisensory 2', 'Unisensory 3', 'Unisensory 4','Unisensory 5','Unisensory 6')
    ylim([-.08, 0.1])
    xlabel('Time (ms)')
    ylabel('Pupil dilation magnitude')
    ax = gca; 
    ax.FontSize = 16; 
    
end

%% this script was first drafted by using data from 4 subjects

% id_name="P29"

% id_name="P32"

% id_name="P49"

% id_name="P50"


%% Up-sample to 1000Hz (not implemented)
%     % Copied from pupil-size-master code (Kret, & Sjak-Shie, 2019)
%     interp_upsamplingFreq = 1000;
%     % Generate the upsampled time vector (seconds):
%     t_upsampled = ...
%         (time(1)/1000 ...
%         :(1/interp_upsamplingFreq)...
%         :time(end)/1000)';
%     diaInterp = interp1(time./1000 ...
%         ,pupilL...
%         ,t_upsampled,'linear'); % this does not interpolate missing data
%     % the resulting vector diaInterp has a frequency of 1000 Hz but
%     % has missing data

%% Filter the data (not implemented)

%    % Low pass filter
%    % Copied from pupil-size-master code (Kret, & Sjak-Shie, 2019)
%    % Calculate the low pass filter specs using the cutoff
%    % frequency [Hz], filter order, and the upsample frequency
%    % specified above:
%    LpFilt_cutoffFreq         = 4;
%    LpFilt_order              = 4;
%    [LpFilt_B,LpFilt_A] ...
%        = butter(LpFilt_order,2*LpFilt_cutoffFreq...
%        /interp_upsamplingFreq );
%    diaInterp_filtered = filtfilt(LpFilt_B...
%        ,LpFilt_A...
%        ,diaInterp);
%     % does not work: error in filtfilt function
%
%     % Multi-step filter by (Kret, & Sjak-Shie, 2019)
%     [valOut,speedFiltData,devFiltData] = rawDataFilter(time, pupilL); #
%     % function from pupil-size-master (Kret, & Sjak-Shie, 2019)
%     % does not work: error in filtfilt function

%% functions needed

function [sRate, msPerS] = etDetermineSampleRate(timeBuffer)
%AMP copied this function from CBCD Luke Mason's database (used in TABLET,
% Face pop-out, probably other Basis tasks)

if isempty(timeBuffer)
    error('Time buffer is empty.')
end

totTimeMs = double(timeBuffer(end, 1) - timeBuffer(1, 1)) / 1000000;
sRate = floor(1 / (totTimeMs / size(timeBuffer, 1)));
msPerS = (totTimeMs / size(timeBuffer, 1)) * 1000;

end

function [dataOut, flags] = etInterpBuffer(lx, ly, rx, ry, timeBuffer, maxMs, dontFilter)

% AMP copied this function from CBCD Luke Mason's database (used in TABLET,
% Face pop-out, probably other Basis tasks)

if ~exist('dontFilter', 'var') || isempty(dontFilter)
    dontFilter = false;
end


flags = false(size(timeBuffer, 1), 4);

% mb = mainBuffer; clear mainBuffer;
tb = timeBuffer; clear timeBuffer;

% get rid of invalid/offscreen data
if ~dontFilter
    %mb = etFilterGazeOnscreen(mb);
end

% get timing data on buffer
[~, msPerS] = etDetermineSampleRate(tb*1000000);
maxSamp = round(maxMs / msPerS);

% % get gaze data
% lx = mb(:, 7);
% ly = mb(:, 8);
% rx = mb(:, 20);
% ry = mb(:, 21);

% find nans
lx_nan = isnan(lx);
ly_nan = isnan(ly);
rx_nan = isnan(rx);
ry_nan = isnan(ry);

% interpolate if there is more than one valid and one invalid sample
t = 1:numel(lx);
if sum(isnan(lx)) > 0 && sum(~isnan(lx)) > 1
    lx_int = interp1(t(~lx_nan), lx(~lx_nan), t, 'linear');
else
    lx_int = lx;
end

if sum(isnan(ly)) > 0 && sum(~isnan(ly)) > 1
    ly_int = interp1(t(~ly_nan), ly(~ly_nan), t, 'linear');
else
    ly_int = ly;
end

if sum(isnan(rx)) > 0 && sum(~isnan(rx)) > 1
    rx_int = interp1(t(~rx_nan), rx(~rx_nan), t, 'linear');
else
    rx_int = rx;
end

if sum(isnan(ry)) > 0 && sum(~isnan(ry)) > 1
    ry_int = interp1(t(~ry_nan), ry(~ry_nan), t, 'linear');
else
    ry_int = ry;
end

% lx_int = interp1q(t(~lx_nan), lx(~lx_nan), t);
% ly_int = interp1q(t(~rx_nan), lx(~rx_nan), t);
% rx_int = interp1q(t(~ly_nan), lx(~ly_nan), t);
% ry_int = interp1q(t(~ry_nan), lx(~ry_nan), t);

% replace interpolated data where the number of missing samples is
% greater than the maximum specified
lx_idx = etInterp_makeIdx(lx, lx_nan, maxSamp);
lx_out = lx;
lx_out(lx_idx) = lx_int(lx_idx);

ly_idx = etInterp_makeIdx(ly, ly_nan, maxSamp);
ly_out = ly;
ly_out(ly_idx) = ly_int(ly_idx);

rx_idx = etInterp_makeIdx(rx, rx_nan, maxSamp);
rx_out = rx;
rx_out(rx_idx) = rx_int(rx_idx);

ry_idx = etInterp_makeIdx(ry, ry_nan, maxSamp);
ry_out = ry;
ry_out(ry_idx) = ry_int(ry_idx);

% store samples that were interpolated in flags output var
flags = any([lx_idx; ly_idx; rx_idx; ry_idx], 1);

% store in buffer
%dataOut = mb;
dataOut = [lx_out, ly_out, rx_out, ry_out];

end

function idx = etInterp_makeIdx(x, xnan, maxSamp)
idx = false(1, length(x));

% find contigous runs of nans, measure length of each
tmp = findcontig(xnan); %, 1);

% if none found, return
if isempty(tmp)
    return
else
    % select those runs that are shorter than the maximum length that
    % we'll interpolate for (default is 150ms)
    tmp = tmp(tmp(:, 3) <= maxSamp, :);
    
    % if none found, return empty
    if isempty(tmp), idx = []; end
    
    % otherwise loop through and mark as valid those segments that are
    % short enough to interpolate
    for e = 1:size(tmp, 1)
        idx(tmp(e, 1):tmp(e, 2)) = true;
    end
end
end

function idx = findcontig(xnan)

idx = zeros(size(xnan,1),3);
a=1;

for ii=1:length(xnan)
    n=1;
    if (ii>1 && xnan(ii-1))
        continue
    end
    if xnan(ii) && (ii+n)<=length(xnan)
        idx(a,1)=ii;
        while xnan(ii+n)
            n=n+1;
            if (ii+n)>length(xnan)
                break
            end
        end
        idx(a,2)=ii+n-1;
        idx(a,3)=n;
        a=a+1;
     end
end

idx = idx(idx(:,1)>0,:);

end
