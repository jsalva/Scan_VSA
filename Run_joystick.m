function Run_VSA(subject_id,run,version,counterbalance,reverse_y_bool)

Screen('Preference', 'SkipSyncTests', 1)

addpath ./util
addpath ./data
KbName('UnifyKeyNames');
debug=false;
global last_refresh;
last_refresh = GetSecs;

filename = ['./data/',subject_id,'_run',num2str(run),'_ver',num2str(version),'_cb',num2str(counterbalance),'_data.txt'];
fileid = fopen(filename,'w');
headerstring = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fileid,headerstring,'Trial','Type','Target/Home(1,0)','Target #','Time(s)','Xpos','Ypos','Xvel','Yvel','Xaccel','Yaccel','Xerror','Yerror','Xprojection','Yprojection');

summary_file = ['./data/',subject_id,'_run',num2str(run),'_ver',num2str(version),'_cb',num2str(counterbalance),'_summary.txt'];
summaryfile = fopen(summary_file,'w');
summaryheader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(summaryfile,summaryheader,'Trial','Type','Target #','Time to Target','Endpoint Error','Est. Avg. Error Mag','Accum. Error x dMagProj', 'Max Error Mag.','Max Error Pos. (% of trajectory)','Feedback','Reached Target');

post_file = ['./data/',subject_id,'_run',num2str(run),'_ver',num2str(version),'_cb',num2str(counterbalance),'_post.txt'];
postfile = fopen(post_file,'w');

event_log = ['./data/',subject_id,'_run',num2str(run),'_ver',num2str(version),'_cb',num2str(counterbalance),'_event.log'];
eventlog = fopen(event_log,'w');
key_format = '%s\t%.4f\n';
notice_format = '%s\n';

%setup experiment paramters
NUMBER_OF_TARGETS = 4;
PERTURBATION_ANGLE = pi()/4;
LENGTH_OF_SEQUENCE = 8;
TARGET_DIST_FROM_CENTER_MM = 100;
TARGET_RADIUS_MM = 5;
WAIT_TIME_POST_TARGET=.25;
%used to setup end point error for positive vs negative feedback
MAX_SUCCESS_THETA = atan2(3*TARGET_RADIUS_MM,TARGET_DIST_FROM_CENTER_MM);

%sound stuff
InitializePsychSound;
[pleasant_data,freq,bits] = wavread('./util/pleasant.wav');
pleasant_data = pleasant_data';
[unpleasant_data,freq,bits] = wavread('./util/unpleasant.wav');
unpleasant_data = 1/16 * [unpleasant_data';unpleasant_data'];
audiohandle=PsychPortAudio('Open',[],[],1);


%experiment block structure and sequences
if ~debug
    %the # of repeats and the trial_duration should result in block
    %durations that are equal across the two runs
    if (counterbalance == 1 && run == 1)||(counterbalance == 2 && run == 2)
        BLOCK_STRUCTURE = {'RD','SD','RD','SD','RD','SD','RD','SD','RD'};
        NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [2,6,2,6,2,6,2,6,2];
        CUE_DURATION = 0.5;
        RANDOM_TRIAL_DURATION = 1.2;
        OTHER_TRIAL_DURATION = 1.2;
        ITI_DURATION = 0.2;
        REST_DURATION = 12;
    elseif (counterbalance == 1 && run == 2)||(counterbalance == 2 && run == 1)
        BLOCK_STRUCTURE = {'RD','RI','RD','RI','RD','RI','RD','RI','RD'};
        NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [2,4,2,4,2,4,2,4,2];
        CUE_DURATION = 0.5;
        RANDOM_TRIAL_DURATION = 1.2;
        OTHER_TRIAL_DURATION = 1.9;
        ITI_DURATION = 0.2;
        REST_DURATION = 12;
    end
    
else
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [1/8 1/8 1/8 1/8 1/8];
    BLOCK_STRUCTURE = {'RD','SD','RD','RI','RD'};
    CUE_DURATION = 0.5;
    OTHER_TRIAL_DURATION = 1.0;
    RANDOM_TRIAL_DURATION = 1.0;
    ITI_DURATION = 0.2;
    REST_DURATION = 1.0;
end
RAND = false;

while (~exist('SEQUENCE','var')) || SEQUENCE(1) == SEQUENCE(end)
    if RAND
        SEQUENCE = sequence_randomizer(LENGTH_OF_SEQUENCE,NUMBER_OF_TARGETS);
    else
        SEQUENCE = [2,4,1,3,4,2,3,1]; %[2 3 1 4 3 2 4 1];% 3 4 2 1];
    end
end

clock = 0;
event_counter = 0;
EVENT_QUEUE = struct();


trial = 1;
toggle = 1;
for block = 1:length(BLOCK_STRUCTURE)

    %Add cue at start of block to the event queue
    event_counter = event_counter + 1;
    EVENT_QUEUE(event_counter).TYPE = 'cue';
    EVENT_QUEUE(event_counter).ONSET = clock;
    clock = clock + CUE_DURATION;
    EVENT_QUEUE(event_counter).TIMEOUT = clock;

    %determine if sequence vs not
    if strcmp(BLOCK_STRUCTURE{block},'SD')
        temp_sequence = repmat(SEQUENCE,1,NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block));
        for i = 1:length(temp_sequence)
            TRIAL(trial).target = temp_sequence(i);
            TRIAL(trial).perturb = false;
            TRIAL(trial).sequence = true;
            TRIAL(trial).type = BLOCK_STRUCTURE{block};

            %add trial to event queue
            event_counter = event_counter + 1;
            EVENT_QUEUE(event_counter).TYPE = 'trial';
            EVENT_QUEUE(event_counter).ONSET = clock;
            clock = clock + OTHER_TRIAL_DURATION;
            EVENT_QUEUE(event_counter).TIMEOUT = clock;

            %add iti to event queue
            event_counter = event_counter + 1;
            EVENT_QUEUE(event_counter).TYPE = 'iti';
            EVENT_QUEUE(event_counter).ONSET = clock;
            clock = clock + ITI_DURATION;
            EVENT_QUEUE(event_counter).TIMEOUT = clock;

            trial = trial + 1;
        end
    %not sequence trial, gen semirandom
    %get temporary sequence
    else
        if RAND
            temp_sequence = sequence_randomizer(LENGTH_OF_SEQUENCE*NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block),NUMBER_OF_TARGETS);
        else
            switch toggle
                case 1
                    temp_sequence = [3,4,1,2,3,1,4,2,3,4,1,2,1,3,2,4,1,2,3,4,2,3,4,1,3,4,2,1,2,4,3,1,4,1,2,3,4,3,1,2,4,1,3,2,3,2,4,1];
                    toggle = 2;
                case 2
                    temp_sequence = [2,3,4,1,4,3,1,2,4,2,1,3,2,3,4,1,2,4,1,3,4,2,3,1,4,1,3,2,4,3,1,2,1,4,3,2,3,4,2,1,4,2,1,3,2,1,3,4];
                    toggle = 3;
                case 3
                    temp_sequence = [1,4,2,3,4,3,1,2,1,2,3,4,4,2,3,1,3,4,2,1,1,4,3,2,2,3,4,1,3,2,1,4,1,4,2,3,1,2,3,4,3,1,4,2,1,3,2,4];
                    toggle = 4;
                case 4
                    temp_sequence = [3,4,2,1,2,4,3,1,4,3,2,1,3,1,4,2,2,4,1,3,2,4,3,1,4,2,3,1,3,4,1,2,1,2,4,3,1,2,4,3,3,4,2,1,1,2,3,4];
                    toggle = 5;
                case 5 
                    temp_sequence = [3,1,4,2,3,4,2,1,2,4,1,3,1,2,4,3,1,2,3,4,1,4,3,2,4,1,2,3,2,3,4,1,2,1,4,3,4,3,1,2,4,3,1,2,4,2,3,1];
                    toggle = 6;
                case 6
                    temp_sequence = [1,2,3,4,1,3,2,4,3,2,1,4,3,4,2,1,4,1,3,2,1,3,2,4,3,1,4,2,4,1,2,3,2,1,3,4,4,1,3,2,4,2,1,3,2,3,1,4];
                    toggle = 7;
                case 7
                    temp_sequence = [3,1,2,4,1,3,4,2,3,1,2,4,3,2,4,1,3,1,2,4,3,1,2,4,2,3,1,4,2,4,3,1,3,4,2,1,4,2,1,3,2,1,4,3,1,3,2,4];
                    toggle = 8;
                case 8
                    temp_sequence = [4,2,3,1,3,1,2,4,1,2,3,4,2,4,1,3,4,1,2,3,1,4,2,3,2,4,1,3,1,3,4,2,2,1,3,4,1,3,2,4,2,1,3,4,4,1,3,2];
                    toggle = 9;
                case 9
                    temp_sequence = [1,2,3,4,3,4,2,1,4,2,3,1,4,3,1,2,3,1,4,2,1,2,3,4,1,3,2,4,2,4,1,3,1,3,4,2,3,2,4,1,1,3,4,2,4,3,1,2];
                    toggle = 10;
                case 10
                    temp_sequence = [1,3,4,2,4,2,3,1,3,4,1,2,1,2,4,3,4,2,1,3,1,2,3,4,2,1,4,3,2,1,3,4,2,4,1,3,4,3,1,2,3,2,4,1,3,1,2,4];
                    toggle = 11;
                case 11
                    temp_sequence = [3,4,1,2,4,1,2,3,3,2,1,4,2,3,1,4,3,2,4,1,4,3,2,1,3,4,2,1,4,2,3,1,2,4,1,3,1,4,2,3,4,3,2,1,2,1,3,4];
                    toggle = 12;
                case 12
                    temp_sequence = [4,2,1,3,1,4,3,2,4,2,1,3,4,3,2,1,2,3,4,1,1,2,4,3,4,3,2,1,4,3,2,1,4,2,1,3,1,3,2,4,2,3,4,1,4,1,3,2];
                    toggle = 13;
                case 13
                    temp_sequence = [4,1,3,2,1,3,4,2,1,2,3,4,1,3,4,2,3,4,2,1,3,1,2,4,1,2,4,3,4,3,2,1,2,3,4,1,4,3,2,1,4,3,1,2,1,4,2,3];
                    toggle = 14;
                case 14
                    temp_sequence = [2,1,4,3,4,2,3,1,2,3,4,1,4,3,2,1,3,4,1,2,1,4,2,3,2,3,1,4,3,2,4,1,2,1,3,4,3,2,4,1,4,2,3,1,3,2,4,1];
                    toggle = 15;
                case 15
                    temp_sequence = [4,2,3,1,3,4,2,1,4,3,2,1,3,1,2,4,3,1,4,2,4,3,2,1,3,2,1,4,1,3,2,4,3,1,4,2,3,4,1,2,3,2,4,1,4,2,1,3];
                    toggle = 16;
                case 16
                    temp_sequence = [1,4,3,2,1,3,4,2,3,1,2,4,2,1,4,3,4,2,3,1,3,2,4,1,4,1,3,2,4,2,1,3,1,4,2,3,4,1,2,3,2,3,4,1,1,4,3,2];
                    toggle = NaN;
            end
        end

        for i = 1:LENGTH_OF_SEQUENCE*NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block)
            TRIAL(trial).target = temp_sequence(i);
            if strcmp(BLOCK_STRUCTURE{block},'RI')
                TRIAL(trial).perturb = true;
            else
                TRIAL(trial).perturb = false;
            end
            TRIAL(trial).sequence = false;
            TRIAL(trial).type = BLOCK_STRUCTURE{block};

            %add trial to event queue
            event_counter = event_counter + 1;
            EVENT_QUEUE(event_counter).TYPE = 'trial';
            EVENT_QUEUE(event_counter).ONSET = clock;
            if strcmp(BLOCK_STRUCTURE{block},'RI')
                clock = clock + OTHER_TRIAL_DURATION;
            else
                clock = clock + RANDOM_TRIAL_DURATION;
            end
            EVENT_QUEUE(event_counter).TIMEOUT = clock;

            %add iti to event queue
            event_counter = event_counter + 1;
            EVENT_QUEUE(event_counter).TYPE = 'iti';
            EVENT_QUEUE(event_counter).ONSET = clock;
            clock = clock + ITI_DURATION;
            EVENT_QUEUE(event_counter).TIMEOUT = clock;

            trial = trial + 1;
        end
    end

    %Add rest at end of block to the event queue
    event_counter = event_counter + 1;
    EVENT_QUEUE(event_counter).TYPE = 'rest';
    EVENT_QUEUE(event_counter).ONSET = clock;
    clock = clock + REST_DURATION;
    EVENT_QUEUE(event_counter).TIMEOUT = clock;

end
event_counter = event_counter + 1;
EVENT_QUEUE(event_counter).TYPE = 'stop';
EVENT_QUEUE(event_counter).ONSET = clock;
EVENT_QUEUE(event_counter).TIMEOUT = clock;
EXPECTED_TTL_COUNT = ceil(clock/2)+1;

TTL_QUEUE = [];
ttl_clock = 0;
for ttl_counter = 1:EXPECTED_TTL_COUNT
    TTL_QUEUE(ttl_counter) = ttl_clock;
    ttl_clock = ttl_clock + 2.0;
end

HideCursor;
global downmap
[keyIsDown sec keyCode] = KbCheckM;
downmap = zeros(size(keyCode));

% Open mainWin on display monitor
screenNum = max(Screen('Screens'));
mainWin = Screen('OpenWindow', screenNum);
mainRect = Screen('Rect', mainWin);

%initialize colors
black = BlackIndex(mainWin);
white = WhiteIndex(mainWin);
grass = [97, 112, 68];

%create main screen
Screen('FillRect', mainWin, grass);

%build cursor screen
cursor_file='./util/golfball.jpg';
cursor_img = imread(cursor_file);
cursor_tex = Screen('MakeTexture', mainWin, cursor_img);
clear cursor_img

%target screen
target_file='./util/golfhole.jpg';
target_img = imread(target_file);
target_tex = Screen('MakeTexture', mainWin, target_img);
clear target_img

%setup stimulus target array
delta_theta = 2*pi()/NUMBER_OF_TARGETS;

%specify screen properties
screen_properties.width_mm = mainRect(3)/4;
screen_properties.height_mm = mainRect(4)/4;
screen_properties.width_res_pix = mainRect(3);
screen_properties.height_res_pix = mainRect(4);
screen_properties.subj_distance_mm = 300;
screen_properties.origin = mainRect/2;
%draw stims

subtends.screen_properties = screen_properties;
%use with rect_subtend_vis_angle
subtends.stim_width_rad = 3.3*2*pi()/360;
subtends.stim_height_rad = 4.8*2*pi()/360;
%use with rect_subtend_distance
subtends.stim_width_mm = TARGET_RADIUS_MM*2;
subtends.stim_height_mm = TARGET_RADIUS_MM*2;

for i = 1:NUMBER_OF_TARGETS
    %% position of each target
    %save cartesian coords (x,y)
    target(i).polar_coords = [i*delta_theta, TARGET_DIST_FROM_CENTER_MM];
    target(i).cartesian_coords = polar2cartesian(screen_properties,i*delta_theta, TARGET_DIST_FROM_CENTER_MM);
    %convert to screen position [0topleft 0topleft x_right y_down] using
    %center as origin
    target(i).screen_position = screen_properties.origin + [0 0 target(i).cartesian_coords(1) -target(i).cartesian_coords(2)];

    subtends.center_on_position = target(i).screen_position;
    target(i).rect = rect_subtend_distance_mm(subtends);
end

%initialize 'home' position (origin in this case)
home.cartesian_coords = [0 0];
home.screen_position = screen_properties.origin;
subtends.center_on_position = home.screen_position;
home.rect = rect_subtend_distance_mm(subtends);

%initialize position object
position.screen_properties=screen_properties;
SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
position.rect = [0 0 0 0];
position.xyt = [0 0 0];
position.polar = [0 0];
position.velxyt = [0 0 0];
position.accelxyt = [0 0 0];
position.error_vec = [0 0];
position.proj_onto_targ = [0 0];
position.error_theta = [0];
position.end_point_error = [0];
position.cue_onset_machine = [0];
position.cue_onset = [0];
position.trial_onset = [0];
position.trial_onset_machine = [0];
position.iti_onset = [0];
position.iti_onset_machine = [0];
position.rest_onset = [0];
position.rest_onset_machine = [0];
position.error_x_dproj_vec= [0];
position.dproj_vec = [0];
position.magerror_vec = [0];
position.magproj_vec = [0];




%ititialize cursor properties
cursor_subtends.screen_properties = screen_properties;
cursor_subtends.stim_width_mm = 5;
cursor_subtends.stim_height_mm = 5;

trial_type = [0];



%setup catch screen
%%%%%%%%%%%%%%%%%%%%
subtends.stim_width_mm = 20;
subtends.stim_height_mm = 20;
subtends.center_on_position=screen_properties.origin + [0 0 0 .2*screen_properties.height_res_pix];
stationary_rect = rect_subtend_distance_mm(subtends);

scale_width = -.2;

TRIGGER_KEY = KbName('=+');
QUIT_KEY = KbName('q');
ttl_counter = 1;

%Create Queue
KbQueueCreate(6)
KbQueueStart()
%%%%%%%%%%%GET QUEUE SET UP%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(eventlog,notice_format,['::Queue started at ',datestr(now),', with machine clock reading ',num2str(GetSecs),'.::']);

while(1)
    
    
    queue_time = GetSecs;
    [pressed firstPress firstRelease lastPress lastRelease] = KbQueueCheck();
    if pressed
        key_ids = find(firstPress);
        first_times = firstPress(key_ids)+ queue_time;
        last_times = lastPress(key_ids) + queue_time;
        key_queue = [0,0];
        queue_counter = 1;
        for key = 1:length(key_ids)
            key_queue(queue_counter,:) = [key_ids(key),first_times(key)];
            queue_counter = queue_counter + 1;

            if first_times(key) ~= last_times(key)
                key_queue(queue_counter,:) = [key_ids(key),last_times(key)];
                queue_counter = queue_counter + 1;
            end
        end
        time_sorted_key_queue = sortrows(key_queue,2);
        for i = 1:length(time_sorted_key_queue(:,1))
            if time_sorted_key_queue(i,1) == TRIGGER_KEY
                position.t0 = GetSecs;%time_sorted_key_queue(i,2);
                fprintf(eventlog,notice_format,['::Experiment started at ',datestr(now),', with machine clock reading ',num2str(position.t0),'.::']);
                fprintf(eventlog,notice_format,['::TTL #',num2str(ttl_counter),' expected: ',num2str(TTL_QUEUE(ttl_counter)),' recorded: 0 ::']);
                ttl_counter = ttl_counter + 1;
                fprintf(eventlog,key_format,KbName(time_sorted_key_queue(i,1)),time_sorted_key_queue(i,2));
                
            else
                fprintf(eventlog,key_format,KbName(time_sorted_key_queue(i,1)),time_sorted_key_queue(i,2));
            end
        end

        if any(key_ids == TRIGGER_KEY)
           break; 
        end
        
        if any(key_ids == QUIT_KEY)
            fprintf(eventlog,notice_format,['::Experimenter quit at ',datestr(now),', with machine clock reading ',num2str(GetSecs),'.::']);
            sca;
            return;
        end
    end
    subtends.center_on_position = screen_properties.origin + [0 0 scale_width*screen_properties.width_res_pix .2*screen_properties.height_res_pix];
    moving_rect = rect_subtend_distance_mm(subtends);
    scale_width = scale_width - .02*scale_width;
    if( scale_width > -.001)
        scale_width = -.2;
    end
    Screen('TextSize',mainWin,100);
    DrawFormattedText(mainWin,'~Quantum Golf~','center',.1*screen_properties.height_res_pix,black);
    Screen('TextSize',mainWin,50);
    
    DrawFormattedText(mainWin,'Goal:','center',.4*screen_properties.height_res_pix,black);
    
    Screen('DrawTexture',mainWin,cursor_tex,[],moving_rect);
    Screen('DrawTexture',mainWin,target_tex,[],stationary_rect);
    Screen('Flip',mainWin);
    
end

freeze = false;
SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
trial = 1;
event_counter = 1;
onset_iteration = true;
within_time_limit = 0;

while(~strcmp(EVENT_QUEUE(event_counter).TYPE,'stop'))
    switch EVENT_QUEUE(event_counter).TYPE
        case 'cue'
            %draw the cue screen depending on whether it's a Perturbed or Not trial
            if TRIAL(trial).perturb
                Screen('FillRect', mainWin, black, home.rect);
            else
                Screen('FillOval',mainWin, black, home.rect);
            end

        case 'trial'
            if ~freeze
                target_num = TRIAL(trial).target;
                %get polar 'theta' for current target, and get current position
                if TRIAL(trial).perturb
                    position.target_theta = target(target_num).polar_coords(1) - PERTURBATION_ANGLE;
                    if reverse_y_bool
                        position = update_pos_data_perturb_reverse_y(position,PERTURBATION_ANGLE);
                    else
                        position = update_pos_data_perturb(position,PERTURBATION_ANGLE);
                    end
                else
                    position.target_theta = target(target_num).polar_coords(1);
                    if reverse_y_bool
                        position = update_pos_data_reverse_y(position);    
                    else
                        position = update_pos_data(position);
                    end
                end
                %convert mouse position to polar
                mag_d = position.polar(size(position.polar,1),2);
                theta_d = position.polar(size(position.polar,1),1);
                %calculate cursor rect
                cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
                cursor_rect = rect_subtend_distance_mm(cursor_subtends);
                %draw cursor and target
                Screen('DrawTexture',mainWin,target_tex,[],target(target_num).rect);
                Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);
                %handle target reached early:
                if mag_d > TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM 
                    within_time_limit = 1;
                    if strcmp(TRIAL(trial).type,'RD')
                        trial_type = 1;
                    elseif strcmp(TRIAL(trial).type,'SD')
                        trial_type = 2;
                    elseif strcmp(TRIAL(trial).type,'RI')
                        trial_type = 3;
                    end
                    
                    time = GetSecs - position.t0;
                    position.time_to_target_from_trial_onset(trial) = time - position.trial_onset(trial);
                    position.end_point_error(trial) = position.error_theta(length(position.error_theta));
                    if abs(position.end_point_error(trial)) < MAX_SUCCESS_THETA
                       PsychPortAudio('FillBuffer',audiohandle,pleasant_data);
                       PsychPortAudio('Start',audiohandle,1,0,1);
                       feedback = 1;
                    else
                       PsychPortAudio('FillBuffer',audiohandle,unpleasant_data);
                       PsychPortAudio('Start',audiohandle,1,0,1);
                       feedback = 0;
                    end

                    positive_error = abs(position.magerror_vec);
                    [max_error_mag,max_index] = nanmax(positive_error);
                    positive_d_projection = [0,abs(diff(position.magproj_vec))];
                    if nansum(positive_d_projection) > 0
                        max_error_pos_as_percent = 100 * nansum(positive_d_projection(1:max_index)) / nansum(positive_d_projection);
                        error_accumulated_over_trajectory = nansum(positive_error .* positive_d_projection );
                        avg_error_mag = error_accumulated_over_trajectory / nansum(positive_d_projection);
                    else
                        max_error_pos_as_percent = 0;
                        error_accumulated_over_trajectory = 0;
                        avg_error_mag = 0;
                    end
                    formatstring = '%d\t%d\t%d\t%d\t%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n';
                    fprintf(fileid,formatstring,[repmat(trial,length(position.xyt(2:end,3)),1),repmat(trial_type,length(position.xyt(2:end,3)),1),(position.xyt(2:end,3)>position.trial_onset(trial)),repmat(target_num,length(position.xyt(2:end,3)),1),position.xyt(2:end,3),position.xyt(2:end,1:2), position.velxyt(2:end,1:2), position.accelxyt(2:end,1:2), position.error_vec(2:end,1:2), position.proj_onto_targ(2:end,1:2)]');
                        
                    summaryformatstring = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n';

                    fprintf(summaryfile,summaryformatstring,[trial,trial_type,target_num,position.time_to_target_from_trial_onset(trial),position.end_point_error(trial),avg_error_mag,error_accumulated_over_trajectory,max_error_mag,max_error_pos_as_percent,feedback,within_time_limit]);

                    position.xyt = [0 0 0];
                    position.velxyt = [0 0 0];
                    position.accelxyt = [0 0 0];
                    position.error_vec = [0 0];
                    position.proj_onto_targ = [0 0];
                    trial_type = [0];
                    position.error_theta = [0];
                    position.rect = [0 0 0 0];
                    position.polar = [0 0];
                    position.error_x_dproj_vec = [0];
                    position.dproj_vec = [0];
                    position.magerror_vec = [0];
                    position.magproj_vec = [0];
                    freeze = true;

                end
            else
                %handle leaving the ball where it ended
                Screen('DrawTexture',mainWin,target_tex,[],target(target_num).rect);
                Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);
            end
        case 'iti'
            %center the ball
            cursor_subtends.center_on_position = screen_properties.origin;
            cursor_rect = rect_subtend_distance_mm(cursor_subtends);
            Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);
        case 'rest'
            Screen('TextSize',mainWin,50);
            DrawFormattedText(mainWin,'+','center','center',black);
    end
    %flip screen and if it's the onset, record and store the onset in machine and relative time
    Screen('Flip',mainWin);
    flip_time = GetSecs - position.t0;
    [pressed firstPress firstRelease lastPress lastRelease] = KbQueueCheck();
    if pressed
        key_ids = find(firstPress);
        first_times = firstPress(key_ids) + flip_time;
        last_times = lastPress(key_ids) + flip_time;
        key_queue = [0,0];
        queue_counter = 1;
        for key = 1:length(key_ids)
            key_queue(queue_counter,:) = [key_ids(key),first_times(key)];
            queue_counter = queue_counter + 1;

            if first_times(key) ~= last_times(key)
                key_queue(queue_counter,:) = [key_ids(key),last_times(key)];
                queue_counter = queue_counter + 1;
            end
        end
        time_sorted_key_queue = sortrows(key_queue,2);
        for i = 1:length(time_sorted_key_queue(:,1))
            if time_sorted_key_queue(i,1) == TRIGGER_KEY
                if ttl_counter <= length(TTL_QUEUE)
                    expected_ttl_arrival_time = num2str(TTL_QUEUE(ttl_counter));
                else
                    expected_ttl_arrival_time = inf;
                end
               
                fprintf(eventlog,notice_format,['::TTL #',num2str(ttl_counter),' expected: ',expected_ttl_arrival_time,' recorded: ',num2str(flip_time),' ::']);
                ttl_counter = ttl_counter + 1;
            end
            fprintf(eventlog,key_format,KbName(time_sorted_key_queue(i,1)),flip_time);
        end

        if any(key_ids == QUIT_KEY)
            fprintf(eventlog,notice_format,['::Experimenter quit at ',datestr(now),', with machine clock reading ',num2str(GetSecs - position.t0),'.::']);
            sca;
            return;
        end
    end
    
    if onset_iteration
        
        fprintf(eventlog,notice_format,['::Onset of ',EVENT_QUEUE(event_counter).TYPE,' expected: ',num2str(EVENT_QUEUE(event_counter).ONSET),' recorded: ',num2str(flip_time),' ::']);
        fprintf(eventlog,key_format,EVENT_QUEUE(event_counter).TYPE,flip_time);

        switch EVENT_QUEUE(event_counter).TYPE
            %use these to make a separate onsets file at the end of the run
            %:)
            case 'cue'
                position.cue_onset(trial) = flip_time;
            case 'trial'
                position.trial_onset(trial) = flip_time;
            case 'iti'
                position.iti_onset(trial) = flip_time;
            case 'rest'
                position.rest_onset(trial) = flip_time;
        end        
        onset_iteration = false;
    end

    %switch events when timeout is reached
    if(flip_time) >= EVENT_QUEUE(event_counter).TIMEOUT
        if strcmp(EVENT_QUEUE(event_counter).TYPE,'trial')
           if ~within_time_limit
                if strcmp(TRIAL(trial).type,'RD')
                    trial_type = 1;
                elseif strcmp(TRIAL(trial).type,'SD')
                    trial_type = 2;
                elseif strcmp(TRIAL(trial).type,'RI')
                    trial_type = 3;
                end

                position.time_to_target_from_trial_onset(trial) = flip_time - position.t0;
                position.end_point_error(trial) = position.error_theta(length(position.error_theta));
                if abs(position.end_point_error(trial)) < MAX_SUCCESS_THETA
                   PsychPortAudio('FillBuffer',audiohandle,pleasant_data);
                   PsychPortAudio('Start',audiohandle,1,0,1);
                   feedback = 1;
                else
                   PsychPortAudio('FillBuffer',audiohandle,unpleasant_data);
                   PsychPortAudio('Start',audiohandle,1,0,1);

                   feedback = 0;
                end

                positive_error = abs(position.magerror_vec);
                [max_error_mag,max_index] = nanmax(positive_error);
                positive_d_projection = [0,abs(diff(position.magproj_vec))];
                if nansum(positive_d_projection) > 0
                    max_error_pos_as_percent = 100 * nansum(positive_d_projection(1:max_index)) / nansum(positive_d_projection);
                    error_accumulated_over_trajectory = nansum(positive_error .* positive_d_projection );
                    avg_error_mag = error_accumulated_over_trajectory / nansum(positive_d_projection);
                else
                    max_error_pos_as_percent = 0;
                    error_accumulated_over_trajectory = 0;
                    avg_error_mag = 0;
                end
                formatstring = '%d\t%d\t%d\t%d\t%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n';
                fprintf(fileid,formatstring,[repmat(trial,length(position.xyt(2:end,3)),1),repmat(trial_type,length(position.xyt(2:end,3)),1),(position.xyt(2:end,3)>position.trial_onset(trial)),repmat(target_num,length(position.xyt(2:end,3)),1),position.xyt(2:end,3),position.xyt(2:end,1:2), position.velxyt(2:end,1:2), position.accelxyt(2:end,1:2), position.error_vec(2:end,1:2), position.proj_onto_targ(2:end,1:2)]');
                    
                summaryformatstring = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n';
                fprintf(summaryfile,summaryformatstring,[trial,trial_type,target_num,position.time_to_target_from_trial_onset(trial),position.end_point_error(trial),avg_error_mag,error_accumulated_over_trajectory,max_error_mag,max_error_pos_as_percent,feedback,within_time_limit]);

                position.xyt = [0 0 0];
                position.velxyt = [0 0 0];
                position.accelxyt = [0 0 0];
                position.error_vec = [0 0];
                position.proj_onto_targ = [0 0];
                trial_type = [0];
                position.error_theta = [0];
                position.rect = [0 0 0 0];
                position.polar = [0 0];
                position.error_x_dproj_vec = [0];
                position.dproj_vec = [0];
                position.magerror_vec = [0];
                position.magproj_vec = [0];
            end

            trial = trial + 1;
            freeze = false;

        end
        event_counter = event_counter + 1;
        onset_iteration = true;
        within_time_limit = 0;
        %re-center cursors
        SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);

    end

end

%%Handle post test
if run == 2

    fprintf(eventlog,notice_format,['::Post-Test started at ',datestr(now),', with machine clock reading ',num2str(position.t0),'.::']);
    STOP = KbName('return');
    string = '';
    line_length = 50;
    line_counter = 1;
    disp_char = '';
    alphabet = ['a':'z','.1234567890'];
    abc123 = {};
    for i = 1:length(alphabet)
        abc123{i} = KbName(alphabet(i));
    end
    fprintf(eventlog,notice_format,['::Post text question started at ',datestr(now),', with machine clock reading ',num2str(position.t0),'.::']);

    while(1)
    
        Screen('Flip',mainWin);
        flip_time = GetSecs - position.t0;
        [pressed firstPress firstRelease lastPress lastRelease] = KbQueueCheck();
        if pressed
            key_ids = find(firstPress);
            first_times = firstPress(key_ids) + flip_time;
            last_times = lastPress(key_ids) + flip_time;
            key_queue = [0,0];
            queue_counter = 1;
            for key = 1:length(key_ids)
                key_queue(queue_counter,:) = [key_ids(key),first_times(key)];
                queue_counter = queue_counter + 1;

                if first_times(key) ~= last_times(key)
                    key_queue(queue_counter,:) = [key_ids(key),last_times(key)];
                    queue_counter = queue_counter + 1;
                end
            end
            time_sorted_key_queue = sortrows(key_queue,2);
            for i = 1:length(time_sorted_key_queue(:,1))
                switch time_sorted_key_queue(i,1) 
                    case KbName('1!')
                        disp_char = '1';
                    case KbName('2@')
                        disp_char = '2';
                    case KbName('3#')
                        disp_char = '3';
                    case KbName('4$')
                        disp_char = '4';
                    case KbName('5%')
                        disp_char = '5';
                    case KbName('6^')
                        disp_char = '6';
                    case KbName('7&')
                        disp_char = '7';
                    case KbName('8*')
                        disp_char = '8';
                    case KbName('9(')
                        disp_char = '9';
                    case KbName('0)')
                        disp_char = '0';
                    case {KbName('space'),KbName('tab')}
                        disp_char = ' ';
                        if length(string)/line_counter > line_length
                            disp_char = [disp_char,'\n'];
                            line_counter = line_counter + 1;
                        end
                    case ~ismac && KbName('backspace')
                        disp_char = '';
                        if ~isempty(string)
                            if length(string)>2 && strcmp(string(end-1:end),'\n')
                                string(end-1:end) = '';
                                line_counter = line_counter - 1;
                            else
                                string(end) = '';
                            end
                        end
                    case KbName('delete')
                        disp_char = '';
                        if ~isempty(string)
                            if length(string)>2 && strcmp(string(end-1:end),'\n')
                                string(end-1:end) = '';
                                line_counter = line_counter - 1;
                            else
                                string(end) = '';
                            end
                        end
                    case KbName(',<')
                        disp_char = ',';
                    case KbName('.>')
                        disp_char = '.';
                    case KbName('''"')
                        disp_char = '''';
                    case KbName('/?')
                        disp_char = '?';
                    case abc123
                        disp_char = KbName(time_sorted_key_queue(i,1));
                    otherwise
                        disp_char = '';
                end

                string = [string, disp_char];
                fprintf(eventlog,key_format,KbName(time_sorted_key_queue(i,1)),flip_time);
            end


            if any(key_ids == STOP)
                break;
            end
        end

        Screen('TextSize',mainWin,25);

        DrawFormattedText(mainWin,string,'center','center',black);

        DrawFormattedText(mainWin,'Did you notice anything interesting about the golf targets?\nPress "Enter" when finished typing.','center',.2*screen_properties.height_res_pix,black);

    end
    postformatstring = '%s\t%s\n';

    type_of_post='Guess about targets:';
    fprintf(postfile,postformatstring,type_of_post,string);




    done_with_post_test = false;
    number = [KbName('1') KbName('2') KbName('3') KbName('4') KbName('1!') KbName('2@') KbName('3#') KbName('4$')];
    counter = 0;
    disp_number = '';
    fprintf(eventlog,notice_format,['::Post numerical sequence started at ',datestr(now),', with machine clock reading ',num2str(position.t0),'.::']);

    while(counter < length(SEQUENCE))

        if done_with_post_test
            break;
        end
        
        Screen('Flip',mainWin);
        flip_time = GetSecs - position.t0;
        [pressed firstPress firstRelease lastPress lastRelease] = KbQueueCheck();
        if pressed
            key_ids = find(firstPress);
            first_times = firstPress(key_ids) + flip_time;
            last_times = lastPress(key_ids) + flip_time;
            key_queue = [0,0];
            queue_counter = 1;
            for key = 1:length(key_ids)
                key_queue(queue_counter,:) = [key_ids(key),first_times(key)];
                queue_counter = queue_counter + 1;

                if first_times(key) ~= last_times(key)
                    key_queue(queue_counter,:) = [key_ids(key),last_times(key)];
                    queue_counter = queue_counter + 1;
                end
            end
            time_sorted_key_queue = sortrows(key_queue,2);
            for i = 1:length(time_sorted_key_queue(:,1))
                switch time_sorted_key_queue(i,1) 
                    case KbName('1!')
                        disp_char = '1';
                        counter = counter + 1;
                    case KbName('2@')
                        disp_char = '2';
                        counter = counter + 1;
                    case KbName('3#')
                        disp_char = '3';
                        counter = counter + 1;
                    case KbName('4$')
                        disp_char = '4';
                        counter = counter + 1;
                    case KbName({'1','2','3','4'})
                        disp_char = KbName(time_sorted_key_queue(i,1));
                        counter = counter + 1;
                    case ~ismac && KbName('backspace')
                        disp_char = '';
                        if ~isempty(disp_number)
                            disp_number(end) = '';
                            counter = counter - 1;
                        end
                    case KbName('delete')
                        disp_char = '';
                        if ~isempty(disp_number)
                            disp_number(end) = '';
                            counter = counter - 1;
                        end
                    otherwise
                        disp_char = '';
                end
                disp_number = [disp_number,disp_char];
            end
        end

        for i = 1:NUMBER_OF_TARGETS
            Screen('TextSize',mainWin,50);
            Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
            DrawFormattedText(mainWin,num2str(i),target(i).rect(1)+.25*(target(i).rect(3)-target(i).rect(1)),target(i).rect(2),[255 0 0]);
        end

        Screen('TextSize',mainWin,25);

        DrawFormattedText(mainWin,disp_number,'center','center',black);

        DrawFormattedText(mainWin,'During some trials\n you were presented with a \nsequence of locations. \n\nEnter the 8 item sequence:','center',.2*screen_properties.height_res_pix,black);


    end
    postformatstring = '%s\t%s\n';

    type_of_post='numbers:';
    disp(disp_number)
    fprintf(postfile,postformatstring,type_of_post,disp_number);

    start_post = GetSecs;
    INSTDUR = 5;
    BETWEEN_STIMULI = true;
    while(GetSecs-start_post < INSTDUR)
        Screen('TextSize',mainWin,25);
        DrawFormattedText(mainWin,'During some trials\n you were presented with a \nsequence of locations. \n\nMove through the 8 item sequence','center',.2*screen_properties.height_res_pix,black);
        for i = 1:NUMBER_OF_TARGETS
            Screen('TextSize',mainWin,50);
            Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
        end
        Screen('Flip',mainWin);
    end


    targ_thetas = [];
    for i = 1:length(target)
        targ_thetas(i) = target(i).polar_coords(1);
    end
    counter = 0;
    position.t0_post = GetSecs;
    next=true;
    draw='home';
    disp_number = '';
    
    fprintf(eventlog,notice_format,['::Post movement (non-forced) sequence started at ',datestr(now),', with machine clock reading ',num2str(position.t0),'.::']);
    while(counter < length(SEQUENCE))

        if done_with_post_test
            break;
        end


        for i = 1:NUMBER_OF_TARGETS
            Screen('TextSize',mainWin,50);
            Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
        end


        if BETWEEN_STIMULI
            SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
        end
        if reverse_y_bool
            position = update_pos_data_reverse_y(position);            
        else
            position = update_pos_data(position);
        end
        %get polar coords of current pos
        mag_d = position.polar(size(position.polar,1),2);
        theta_d = position.polar(size(position.polar,1),1);

        %draw cursor centered on the current position
        cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
        cursor_rect = rect_subtend_distance_mm(cursor_subtends);
        %Screen('FillOval',mainWin,white,cursor_rect);
        Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);

        switch draw
            case 'home'
                if next
                    BETWEEN_STIMULI=true;
                    position.home_onset_post(counter+1) = GetSecs - position.t0_post;
                    next = false;
                end

                if GetSecs - position.t0_post - position.home_onset_post(counter+1) > WAIT_TIME_POST_TARGET
                    next = true;
                    BETWEEN_STIMULI=false;
                    draw = 'target';
                end

            case 'target'
                if mag_d > TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM
                    draw = 'home';
                    counter = counter+1;
                    if theta_d<2*pi()/(2*NUMBER_OF_TARGETS)
                        disp_number(counter)=num2str(NUMBER_OF_TARGETS);
                    else
                        disp_number(counter)=num2str(find(abs(targ_thetas-theta_d)==min(abs(targ_thetas-theta_d))));
                    end
                end

        end


        Screen('Flip',mainWin);
        flip_time = GetSecs - position.t0;
        [pressed firstPress firstRelease lastPress lastRelease] = KbQueueCheck();
        if pressed
            key_ids = find(firstPress);
            first_times = firstPress(key_ids) + flip_time;
            last_times = lastPress(key_ids) + flip_time;
            key_queue = [0,0];
            queue_counter = 1;
            for key = 1:length(key_ids)
                key_queue(queue_counter,:) = [key_ids(key),first_times(key)];
                queue_counter = queue_counter + 1;

                if first_times(key) ~= last_times(key)
                    key_queue(queue_counter,:) = [key_ids(key),last_times(key)];
                    queue_counter = queue_counter + 1;
                end
            end
            time_sorted_key_queue = sortrows(key_queue,2);
            for i = 1:length(time_sorted_key_queue(:,1))
                if time_sorted_key_queue(i,1) == TRIGGER_KEY
                    if ttl_counter <= length(TTL_QUEUE)
                        expected_ttl_arrival_time = num2str(TTL_QUEUE(ttl_counter));
                    else
                        expected_ttl_arrival_time = inf;
                    end

                    fprintf(eventlog,notice_format,['::TTL #',num2str(ttl_counter),' expected: ',expected_ttl_arrival_time,' recorded: ',num2str(flip_time),' ::']);
                    ttl_counter = ttl_counter + 1;
                end
                fprintf(eventlog,key_format,KbName(time_sorted_key_queue(i,1)),flip_time);
            end

            if any(key_ids == QUIT_KEY)
                fprintf(eventlog,notice_format,['::Experimenter quit at ',datestr(now),', with machine clock reading ',num2str(GetSecs - position.t0),'.::']);
                sca;
                return;
            end
        end
    end
    disp(disp_number)
    type_of_post='movements:';
    fprintf(postfile,postformatstring,type_of_post,disp_number);

    start_post = GetSecs;
    INSTDUR = 5;
    

    while(GetSecs-start_post < INSTDUR)
        Screen('TextSize',mainWin,25);
        DrawFormattedText(mainWin,'During some trials\n you were presented with a \nsequence of locations. \n\nMove through the 8 item sequence','center',.2*screen_properties.height_res_pix,black);
        for i = 1:NUMBER_OF_TARGETS
            Screen('TextSize',mainWin,50);
            Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
        end
        Screen('Flip',mainWin);
    end

    %Correct sequence needed
    counter = 0;
    disp_tmp=[];
    wait = true;
    disp_number='';
    draw='home';
    next=true;
    fprintf(eventlog,notice_format,['::Post movement (forced) sequence started at ',datestr(now),', with machine clock reading ',num2str(position.t0),'.::']);

    while(counter < length(SEQUENCE))

        if done_with_post_test
            break;
        end


        for i = 1:NUMBER_OF_TARGETS
            Screen('TextSize',mainWin,50);
            Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
        end


        if BETWEEN_STIMULI
            SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
        end

        if reverse_y_bool
            position = update_pos_data_reverse_y(position);            
        else
            position = update_pos_data(position);
        end

        %get polar coords of current pos
        mag_d = position.polar(size(position.polar,1),2);
        theta_d = position.polar(size(position.polar,1),1);

        %draw cursor centered on the current position
        cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
        cursor_rect = rect_subtend_distance_mm(cursor_subtends);
        %Screen('FillOval',mainWin,white,cursor_rect);
        Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);

        switch draw
            case 'home'
                if next
                    BETWEEN_STIMULI=true;
                    position.home_onset_post(counter+1) = GetSecs - position.t0_post;
                    next = false;
                end

                if GetSecs - position.t0_post - position.home_onset_post(counter+1) > WAIT_TIME_POST_TARGET
                    next = true;
                    BETWEEN_STIMULI=false;
                    draw = 'target';
                end

            case 'target'
                if mag_d > TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM  && wait


                    if theta_d<2*pi()/(2*NUMBER_OF_TARGETS)
                        disp_tmp(length(disp_tmp)+1)=NUMBER_OF_TARGETS;
                    else
                        disp_tmp(length(disp_tmp)+1)=find(abs(targ_thetas-theta_d)==min(abs(targ_thetas-theta_d)));
                    end

                    disp_number(length(disp_number)+1) = num2str(disp_tmp(length(disp_tmp)));
                    if disp_tmp(length(disp_tmp)) == SEQUENCE(counter+1)
                        disp_tmp=[];
                        draw = 'home';
                        disp_number(length(disp_number)+1) = '-';
                        counter = counter+1; 
                    end
                    wait=false;
                elseif mag_d < TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM
                    wait = true;
                end

        end


        Screen('Flip',mainWin);
        flip_time = GetSecs - position.t0;
        [pressed firstPress firstRelease lastPress lastRelease] = KbQueueCheck();
        if pressed
            key_ids = find(firstPress);
            first_times = firstPress(key_ids) + flip_time;
            last_times = lastPress(key_ids) + flip_time;
            key_queue = [0,0];
            queue_counter = 1;
            for key = 1:length(key_ids)
                key_queue(queue_counter,:) = [key_ids(key),first_times(key)];
                queue_counter = queue_counter + 1;

                if first_times(key) ~= last_times(key)
                    key_queue(queue_counter,:) = [key_ids(key),last_times(key)];
                    queue_counter = queue_counter + 1;
                end
            end
            time_sorted_key_queue = sortrows(key_queue,2);
            for i = 1:length(time_sorted_key_queue(:,1))
                if time_sorted_key_queue(i,1) == TRIGGER_KEY
                    if ttl_counter <= length(TTL_QUEUE)
                        expected_ttl_arrival_time = num2str(TTL_QUEUE(ttl_counter));
                    else
                        expected_ttl_arrival_time = inf;
                    end

                    fprintf(eventlog,notice_format,['::TTL #',num2str(ttl_counter),' expected: ',expected_ttl_arrival_time,' recorded: ',num2str(flip_time),' ::']);
                    ttl_counter = ttl_counter + 1;
                end
                fprintf(eventlog,key_format,KbName(time_sorted_key_queue(i,1)),flip_time);
            end

            if any(key_ids == QUIT_KEY)
                fprintf(eventlog,notice_format,['::Experimenter quit at ',datestr(now),', with machine clock reading ',num2str(GetSecs - position.t0),'.::']);
                sca;
                return;
            end
        end
    end
    disp(disp_number)
    type_of_post='forced crct:';
    fprintf(postfile,postformatstring,type_of_post,disp_number);
end




sca
%save('./util/cb.mat','counterbal');




end
