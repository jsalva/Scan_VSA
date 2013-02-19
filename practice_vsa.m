function practice_vsa(subject_id,run)

Screen('Preference', 'SkipSyncTests', 1)
if run == 2
    reverse_y_bool=true;
else
    reverse_y_bool=false;
end
addpath ./util
addpath ./data
KbName('UnifyKeyNames');
debug=false;
global last_refresh;
last_refresh = GetSecs;

filename = ['./data/',subject_id,'_pratice',num2str(run),'_data.txt'];
fileid = fopen(filename,'w');
headerstring = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fileid,headerstring,'Trial','Type','Target/Home(1,0)','Target #','Time(s)','Xpos','Ypos','Xvel','Yvel','Xaccel','Yaccel','Xerror','Yerror','Xprojection','Yprojection');

summary_file = ['./data/',subject_id,'_pratice',num2str(run),'_summary.txt'];
summaryfile = fopen(summary_file,'w');
summaryheader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(summaryfile,summaryheader,'Trial','Type','Onset','Target #','Time to Target','Endpoint Error','Est. Avg. Error Mag','Accum. Error x dMagProj', 'Max Error Mag.','Max Error Pos. (% of trajectory)','Feedback','Target (1)/Return (0)');

event_log = ['./data/',subject_id,'pratice',num2str(run),'_event.log'];
eventlog = fopen(event_log,'w');
key_format = '%s\t%.4f\n';
notice_format = '%s\n';

onsets_file = ['./data/',subject_id,'pratice',num2str(run),'_onsets.txt'];
onsetsfile = fopen(onsets_file,'w');
onsetsheader_format = '%s\t%s\t%s\n';
onsets_format = '%d\t%s\t%0.10f\n';
fprintf(onsetsfile,onsetsheader_format,'Event Number','Event Type','Onset');

%setup experiment paramters
NUMBER_OF_TARGETS = 4;
PERTURBATION_ANGLE = pi()/6;
LENGTH_OF_SEQUENCE = 8;
global TARGET_DIST_FROM_CENTER_MM
global JOYSTICK_MAGNITUDE_SCALING

TARGET_DIST_FROM_CENTER_MM = 100;
JOYSTICK_MAGNITUDE_SCALING = 1.3;
RETURN_TARGET_RADIUS_MM = 20;
TARGET_RADIUS_MM = 5;
WAIT_TIME_POST_TARGET=.25;
%used to setup end point error for positive vs negative feedback
MAX_SUCCESS_THETA = atan2(3*TARGET_RADIUS_MM,TARGET_DIST_FROM_CENTER_MM);
TR_DURATION = 2.0;

%sound stuff
InitializePsychSound;
[pleasant_data,freq,bits] = wavread('./util/pleasant.wav');
pleasant_data = pleasant_data';
[unpleasant_data,freq,bits] = wavread('./util/unpleasant.wav');
unpleasant_data = 1/16 * [unpleasant_data';unpleasant_data'];
audiohandle=PsychPortAudio('Open',[],[],1);
if run ==2
    BLOCK_STRUCTURE = {'RD','RD','RD','RD'};
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [6,6,6,6];

    CUE_DURATION = 0.5;

    TRIAL_DURATION = inf;
    TRIAL_TIME_ON_TARGET = 0.1;

    RETURN_DURATION = inf;
    RETURN_TIME_ON_TARGET = 0.1;

    REST_DURATION = 12;
    RAND = true;

else
    BLOCK_STRUCTURE = {'RD','RI'};
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [2,2];

    CUE_DURATION = 0.5;

    TRIAL_DURATION = inf;
    TRIAL_TIME_ON_TARGET = 0.1;

    RETURN_DURATION = inf;
    RETURN_TIME_ON_TARGET = 0.1;

    REST_DURATION = 5;
    
    RAND = false;
end

while (~exist('SEQUENCE','var')) || SEQUENCE(1) == SEQUENCE(end)
    if RAND
        SEQUENCE = sequence_randomizer(LENGTH_OF_SEQUENCE,NUMBER_OF_TARGETS);
    else
        SEQUENCE = [2,4,1,3,4,2,3,1]; %[2 3 1 4 3 2 4 1];% 3 4 2 1];
    end
end

event_counter = 0;
EVENT_QUEUE = struct();


trial = 1;
toggle = 1;
for block = 1:length(BLOCK_STRUCTURE)

    %Add cue at start of block to the event queue
    event_counter = event_counter + 1;
    EVENT_QUEUE(event_counter).TYPE = 'cue';
    EVENT_QUEUE(event_counter).DURATION = CUE_DURATION;

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
            EVENT_QUEUE(event_counter).DURATION = TRIAL_DURATION;
            EVENT_QUEUE(event_counter).TIME_ON_TARGET = TRIAL_TIME_ON_TARGET;
            

            %add return to event queue
            event_counter = event_counter + 1;
            EVENT_QUEUE(event_counter).TYPE = 'return';
            EVENT_QUEUE(event_counter).DURATION = RETURN_DURATION;
            EVENT_QUEUE(event_counter).TIME_ON_TARGET = RETURN_TIME_ON_TARGET;


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
            EVENT_QUEUE(event_counter).DURATION = TRIAL_DURATION;
            EVENT_QUEUE(event_counter).TIME_ON_TARGET = TRIAL_TIME_ON_TARGET;
            

            %add return to event queue
            event_counter = event_counter + 1;
            EVENT_QUEUE(event_counter).TYPE = 'return';
            EVENT_QUEUE(event_counter).DURATION = RETURN_DURATION;
            EVENT_QUEUE(event_counter).TIME_ON_TARGET = RETURN_TIME_ON_TARGET;
            
            trial = trial + 1;
        end
    end

    %Add rest at end of block to the event queue
    event_counter = event_counter + 1;
    EVENT_QUEUE(event_counter).TYPE = 'rest';
    EVENT_QUEUE(event_counter).DURATION = REST_DURATION;

end
event_counter = event_counter + 1;
EVENT_QUEUE(event_counter).TYPE = 'stop';

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
cursor_subtends.center_on_position = home.screen_position;
cursor_home = rect_subtend_distance_mm(cursor_subtends);
cursor_subtends.center_on_position = target(1).screen_position;
cursor_on_target = rect_subtend_distance_mm(cursor_subtends);
off_target = polar2cartesian(screen_properties,1.25*delta_theta, TARGET_DIST_FROM_CENTER_MM);
disp(off_target)
cursor_subtends.center_on_position = screen_properties.origin + [0,0,off_target(1),-off_target(2)];
cursor_off_target = rect_subtend_distance_mm(cursor_subtends);
trial_type = [0];



%setup catch screen
%%%%%%%%%%%%%%%%%%%%
subtends.stim_width_mm = 20;
subtends.stim_height_mm = 20;
subtends.center_on_position=screen_properties.origin + [0 0 0 .2*screen_properties.height_res_pix];
stationary_rect = rect_subtend_distance_mm(subtends);

scale_width = -.2;

TRIGGER_KEY = KbName('space');
QUIT_KEY = KbName('q');
ttl_counter = 1;

%Create Queue
%Currently only works for mac, and is done this way to avoid
%improperly detecting the scanner button box as the queue to open. You will
%have to change this if you're not using a mac.
keyboard = GetKeyboardIndices('Apple Internal Keyboard / Trackpad');
KbQueueCreate(keyboard)
KbQueueStart()
%%%%%%%%%%%GET QUEUE SET UP%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(eventlog,notice_format,['::Queue started at ',datestr(now),', with machine clock reading ',num2str(GetSecs),'.::']);
if run == 1
    t0=GetSecs;
    screen = 1;
    while(1)
        flip_time = GetSecs - t0;
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

            if any(key_ids == QUIT_KEY)
                fprintf(eventlog,notice_format,['::Experimenter quit at ',datestr(now),', with machine clock reading ',num2str(GetSecs - position.t0),'.::']);
                sca;
                return;
            elseif any(key_ids == TRIGGER_KEY)
                if screen == 3
                    PsychPortAudio('FillBuffer',audiohandle,pleasant_data);
                    PsychPortAudio('Start',audiohandle,1,0,1);
                elseif screen == 4
                    PsychPortAudio('FillBuffer',audiohandle,unpleasant_data);
                    PsychPortAudio('Start',audiohandle,1,0,1);
                elseif screen == 7
                    break;
                end
                screen = screen + 1;
            end
        end
        switch screen
            case 1
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'Your goal is to use the joystick to move the golf ball into the target.\nTargets look like the picture above.','center',.2*screen_properties.height_res_pix,black);

                Screen('DrawTexture',mainWin,cursor_tex,[],cursor_home);
                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);
                Screen('DrawTexture',mainWin,target_tex,[],target(1).rect);
            case 2
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'Targets will appear in one of these four locations.\nYour goal is to move the ball onto the target and hold it there until the next target appears.\nGet the ball to the target as quickly as you can, in the most direct, straigt-line path.','center',.2*screen_properties.height_res_pix,black);       

                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);
                for i = 1:NUMBER_OF_TARGETS
                    Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
                end
            case 3
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'If you make it to the target accurately, you will hear a pleasant sound.','center',.2*screen_properties.height_res_pix,black);       
                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);
                Screen('DrawTexture',mainWin,target_tex,[],target(1).rect);
                Screen('DrawTexture',mainWin,cursor_tex,[],cursor_on_target);

            case 4
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'If you fail to make it to the target accurately, you will hear an annoying noise.','center',.2*screen_properties.height_res_pix,black);       
                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);
                Screen('DrawTexture',mainWin,target_tex,[],target(1).rect);
                Screen('DrawTexture',mainWin,cursor_tex,[],cursor_off_target);
            case 5
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'Every other target will be at the center of the screen.','center',.2*screen_properties.height_res_pix,black);       
                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);
                Screen('DrawTexture',mainWin,target_tex,[],home.rect);
            case 6
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'Most of the time, the joystick will move the way you expect.\nHowever, sometimes it will not. When this happens, you will see a black square.\nAs best as you can, continue to move the ball to the target\nas quickly as you can, in the most direct, straigt-line path.','center',.2*screen_properties.height_res_pix,black);       
                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);
                Screen('FillRect', mainWin, black, home.rect);
            case 7
                Screen('TextSize',mainWin,25);
                DrawFormattedText(mainWin,'Get ready to practice for yourself!','center',.2*screen_properties.height_res_pix,black);       
                Screen('TextSize',mainWin,15);
                DrawFormattedText(mainWin,'When you are ready to move on, press the spacebar','center',.8*screen_properties.height_res_pix,black);

        end
        Screen('Flip',mainWin);
    end
end


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
                fprintf(eventlog,notice_format,['::TTL #',num2str(ttl_counter),' expected (time-based): ',num2str(round((GetSecs-position.t0)/TR_DURATION)),' recorded: 0 ::']);
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

%freeze = false;
SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
trial = 1;
event_counter = 1;
onset_iteration = true;
trial_done = false;
iti_done = false;
first_trial_iteration = true;
first_return_iteration = true;
next_event = false;
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
            
            target_num = TRIAL(trial).target;
            
            if TRIAL(trial).perturb
                position.target_theta = target(target_num).polar_coords(1) - PERTURBATION_ANGLE;
                position = update_absolute_pos_data_perturb(position, reverse_y_bool, PERTURBATION_ANGLE);
            else
                position.target_theta = target(target_num).polar_coords(1);
                position = update_absolute_pos_data(position,reverse_y_bool);
            end
            
            mag_d = position.polar(size(position.polar,1),2);
            theta_d = position.polar(size(position.polar,1),1);
            
            cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
            cursor_rect = rect_subtend_distance_mm(cursor_subtends);
            
            Screen('DrawTexture',mainWin,target_tex,[],target(target_num).rect);
            Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);
            
            if mag_d > TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM 
                
                if first_trial_iteration
                   target_onset = GetSecs - position.t0;
                   endpoint = length(position.error_theta);
                   nth_target_onset = target_onset;
                   first_trial_iteration = false;
                end
                
                time_on_target = (GetSecs - position.t0) - nth_target_onset; 
                
                    if time_on_target >= EVENT_QUEUE(event_counter).TIME_ON_TARGET
                        is_target = 1;
                        if strcmp(TRIAL(trial).type,'RD')
                            trial_type = 1;
                        elseif strcmp(TRIAL(trial).type,'SD')
                            trial_type = 2;
                        elseif strcmp(TRIAL(trial).type,'RI')
                            trial_type = 3;
                        end

                        time = GetSecs - position.t0;
                        
                        position.time_to_target_from_trial_onset(trial) = target_onset - position.trial_onset(trial);
                        position.time_to_target_from_trial_onset_final(trial) = time - position.trial_onset(trial);
                        
                        position.end_point_error(trial) = position.error_theta(endpoint);
                        position.end_point_error_final(trial) = position.error_theta(length(position.error_theta));
                        
                        if abs(position.end_point_error_final(trial)) < MAX_SUCCESS_THETA
                           PsychPortAudio('FillBuffer',audiohandle,pleasant_data);
                           PsychPortAudio('Start',audiohandle,1,0,1);
                           feedback = 1;
                        else
                           PsychPortAudio('FillBuffer',audiohandle,unpleasant_data);
                           PsychPortAudio('Start',audiohandle,1,0,1);
                           feedback = 0;
                        end

                        positive_error = abs(position.magerror_vec(1:endpoint));
                        [max_error_mag,max_index] = nanmax(positive_error);
                        
                        positive_d_projection = [0,abs(diff(position.magproj_vec(1:endpoint)))];
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

                        summaryformatstring = '%d\t%d\t%0.4f\t%d\t%0.4f\t%f\t%d\t%d\t%d\t%0.4f\t%d\t%d\n';

                        fprintf(summaryfile,summaryformatstring,[trial,trial_type,position.trial_onset(trial),target_num,position.time_to_target_from_trial_onset(trial),position.end_point_error(trial),avg_error_mag,error_accumulated_over_trajectory,max_error_mag,max_error_pos_as_percent,feedback,is_target]);

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
                        next_event = true;
                        first_trial_iteration = true;

                    end
            
            else
                nth_target_onset = GetSecs - position.t0;
                time_on_target = 0;
            end
    
        case 'return'
            
            target_num = TRIAL(trial).target;
            
            if TRIAL(trial).perturb
                position.target_theta = target(target_num).polar_coords(1) - PERTURBATION_ANGLE;
                position = update_absolute_pos_data_perturb(position, reverse_y_bool, PERTURBATION_ANGLE);
            else
                position.target_theta = target(target_num).polar_coords(1);
                position = update_absolute_pos_data(position,reverse_y_bool);
            end
            
            mag_d = position.polar(size(position.polar,1),2);
            theta_d = position.polar(size(position.polar,1),1);
            
            cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
            cursor_rect = rect_subtend_distance_mm(cursor_subtends);
            
            Screen('DrawTexture',mainWin,target_tex,[],home.rect);
            Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);
            
            if mag_d < RETURN_TARGET_RADIUS_MM
                
                if first_return_iteration
                   return_onset = GetSecs - position.t0;
                   endpoint = length(position.error_theta);
                   nth_return_onset = return_onset;
                   first_return_iteration = false;
                end
                
                time_on_target = (GetSecs - position.t0) - nth_return_onset; 
                
                    if time_on_target >= EVENT_QUEUE(event_counter).TIME_ON_TARGET
                        is_target = 0;
                        trial_type = 0;
                        time = GetSecs - position.t0;
                        
                        time_to_target_from_return_onset = return_onset - position.return_onset(trial);
                        time_to_target_from_return_onset_final = time - position.return_onset(trial);
                        
                        %error is meaningless for return
                        end_point_error = 0;
                        
                        PsychPortAudio('FillBuffer',audiohandle,pleasant_data);
                        PsychPortAudio('Start',audiohandle,1,0,1);
                        feedback = 1;

                        positive_error = abs(position.magerror_vec(1:endpoint));
                        [max_error_mag,max_index] = nanmax(positive_error);
                        
                        positive_d_projection = [0,abs(diff(position.magproj_vec(1:endpoint)))];
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

                        summaryformatstring = '%d\t%d\t%0.4f\t%d\t%0.4f\t%f\t%d\t%d\t%d\t%0.4f\t%d\t%d\n';

                        fprintf(summaryfile,summaryformatstring,[trial,trial_type,position.return_onset(trial),target_num,time_to_target_from_return_onset,end_point_error,avg_error_mag,error_accumulated_over_trajectory,max_error_mag,max_error_pos_as_percent,feedback,is_target]);

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
                        next_event = true;
                        first_return_iteration = true;

                    end
            
            else
                nth_return_onset = GetSecs - position.t0;
                time_on_target = 0;
            end            
            
            
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
                fprintf(eventlog,notice_format,['::TTL #',num2str(ttl_counter),' expected (time-based): ',num2str(round(flip_time/TR_DURATION)),' recorded: ',num2str(flip_time),' ::']);
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
        last_onset = flip_time;
        fprintf(eventlog,notice_format,['::Onset of ',EVENT_QUEUE(event_counter).TYPE,' expected: [N/A SELF PACED] ',' recorded: ',num2str(flip_time),' ::']);
        fprintf(eventlog,key_format,EVENT_QUEUE(event_counter).TYPE,flip_time);

        switch EVENT_QUEUE(event_counter).TYPE
            case 'cue'
                position.cue_onset(trial) = flip_time;
            case 'trial'
                position.trial_onset(trial) = flip_time;
            case 'return'
                position.return_onset(trial) = flip_time;
            case 'rest'
                position.rest_onset(trial) = flip_time;
        end        
        
        fprintf(onsetsfile,onsets_format,event_counter,EVENT_QUEUE(event_counter).TYPE,flip_time);
        onset_iteration = false;
    end

    %switch events when duration is reached (remember it's inf for trial
    %and return trials)
    if next_event || (flip_time) - last_onset >= EVENT_QUEUE(event_counter).DURATION
        
        if strcmp(EVENT_QUEUE(event_counter).TYPE,'return')
            trial = trial + 1;
        end
        onset_iteration = true;
        event_counter = event_counter + 1;
        next_event = false;
        
    end

end

sca
%save('./util/cb.mat','counterbal');




end
