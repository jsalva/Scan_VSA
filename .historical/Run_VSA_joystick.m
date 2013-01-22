function Run_VSA_joystick(subject_id)
Screen('Preference', 'SkipSyncTests', 1)
addpath ./util
addpath ./data
KbName('UnifyKeyNames');
debug=false;
Gamepad;
filename = ['./data/',subject_id,'_data.txt'];
fileid = fopen(filename,'w');
headerstring = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fileid,headerstring,'Trial','Type','Target/Home(1,0)','Target #','Time(s)','Xpos','Ypos','Xvel','Yvel','Xaccel','Yaccel','Xerror','Yerror','Xprojection','Yprojection');
                

summary_file = ['./data/',subject_id,'_summary.txt'];
summaryfile = fopen(summary_file,'w');
summaryheader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(summaryfile,summaryheader,'Trial','Type','Target #','Time to Target','Endpoint Error','Est. Avg. Error Mag','Accum. Error x dMagProj', 'Max Error Mag.','Max Error Pos. (% of trajectory)','Feedback');


post_file = ['./data/',subject_id,'_post.txt'];
postfile = fopen(post_file,'w');


%setup experiment paramters
NUMBER_OF_TARGETS = 4;
PERTURBATION_ANGLE = pi()/4;
LENGTH_OF_SEQUENCE = 8;
TARGET_DIST_FROM_CENTER_MM = 100;
TARGET_RADIUS_MM = 5;
WAIT_TIME_POST_TARGET=.25;
%used to setup end point error for positive vs negative feedback
MAX_SUCCESS_THETA = atan2(3*TARGET_RADIUS_MM,TARGET_DIST_FROM_CENTER_MM);
%cb = load('./util/cb.mat');
%counterbal = cb.counterbal;

%sound stuff
InitializePsychSound;
[pleasant_data,freq,bits] = wavread('./util/pleasant.wav');
pleasant_data = pleasant_data';
[unpleasant_data,freq,bits] = wavread('./util/unpleasant.wav');
unpleasant_data = 1/16 * [unpleasant_data';unpleasant_data'];
audiohandle=PsychPortAudio('Open',[],[],1);


%PsychPortAudio('Volume',audiohandle,.05);

if ~debug
    %if counterbal == 1
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [5 20 5 20 5];
    BLOCK_STRUCTURE = {'RD','SD','RD','RI','RD'};
        %counterbal = 2;
    %elseif counterbal == 2 
    %    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [5 20 5 20 5];
    %    BLOCK_STRUCTURE = {'RD','RI','RD','SD','RD'};
    %    counterbal = 1;
    %end
else
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [1/8 1/8 1/8 1/8 1/8];
    BLOCK_STRUCTURE = {'RD','SD','RD','RI','RD'};
end
RAND = false;

while (~exist('SEQUENCE','var')) || SEQUENCE(1) == SEQUENCE(end)
    if RAND
        SEQUENCE = sequence_randomizer(LENGTH_OF_SEQUENCE,NUMBER_OF_TARGETS);
    else
        SEQUENCE = [2 4 1 3 4 2 3 1]; %[2 3 1 4 3 2 4 1];% 3 4 2 1];
    end
end
trial = 1;
toggle = 1;
for block = 1:length(BLOCK_STRUCTURE)
    %determine if sequence vs not
    if strcmp(BLOCK_STRUCTURE{block},'SD')
	disp(NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block))
        temp_sequence = repmat(SEQUENCE,1,NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block));
        if debug
            TRIAL(trial).target = SEQUENCE(1);
            TRIAL(trial).perturb = false;
            TRIAL(trial).sequence = true;
            TRIAL(trial).type = BLOCK_STRUCTURE{block};
            trial = trial + 1;
        end
        for i = 1:length(temp_sequence)
            TRIAL(trial).target = temp_sequence(i);
            TRIAL(trial).perturb = false;
            TRIAL(trial).sequence = true;
            TRIAL(trial).type = BLOCK_STRUCTURE{block};
            trial = trial + 1;
        end
    %not sequence trial, gen semirandom
    %get temporary sequence
    else
        if RAND
            temp_sequence = sequence_randomizer(LENGTH_OF_SEQUENCE*NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block),NUMBER_OF_TARGETS);
        else
            switch toggle
                %these particular sequences were generated with 8*20 items
                %with 4 targets. Since  as of 5/14/12, there are blocks with 8*10 (1/2), and 8*5 (1/4) the
                %length of the following sequences, counterbalanced across
                %conditions, only the appropriate fraction of the sequences
                %will be used in those blocks
                case 1
                    temp_sequence = [2,3,4,1,2,1,4,3,1,4,2,3,2,1,3,4,3,1,2,4,1,3,2,4,2,3,1,4,3,1,2,4,1,3,2,4,3,4,2,1,3,2,1,4,2,4,3,1,4,2,1,3,2,4,1,3,1,2,3,4,2,1,3,4,2,3,1,4,3,4,1,2,3,4,2,1,4,2,3,1,3,1,2,4,1,2,4,3,2,1,4,3,4,3,1,2,4,1,2,3,1,2,4,3,2,4,3,1,4,3,1,2,1,3,4,2,1,2,3,4,1,4,2,3,1,2,4,3,4,3,2,1,3,4,2,1,3,1,2,4,3,1,2,4,2,1,4,3,1,4,3,2,1,4,2,3,1,4,3,2];
                    toggle = 2;
                case 2
                    temp_sequence = [4,1,2,3,2,4,3,1,3,4,1,2,1,2,4,3,4,3,2,1,4,1,3,2,3,2,4,1,3,2,4,1,3,2,4,1,4,3,2,1,2,4,3,1,4,1,3,2,1,4,2,3,1,3,2,4,3,4,1,2,1,3,4,2,4,3,1,2,1,4,2,3,4,1,2,3,1,3,2,4,2,1,4,3,1,2,4,3,1,2,4,3,1,2,4,3,2,1,3,4,2,4,1,3,4,3,1,2,3,1,4,2,1,3,4,2,1,3,4,2,3,1,2,4,3,2,1,4,1,2,4,3,2,3,4,1,4,1,3,2,3,2,4,1,2,1,3,4,3,4,1,2,1,4,2,3,2,1,3,4];
                    toggle = 3;
                case 3
                    temp_sequence = [1,4,2,3,4,1,3,2,3,2,4,1,2,1,4,3,2,1,3,4,2,4,1,3,1,2,3,4,2,4,1,3,4,2,3,1,3,1,4,2,3,2,1,4,3,4,1,2,3,4,1,2,3,1,2,4,3,1,2,4,2,4,3,1,4,2,1,3,1,2,4,3,2,3,1,4,2,1,4,3,4,2,3,1,3,1,2,4,2,4,1,3,4,3,2,1,4,3,2,1,2,4,1,3,1,3,4,2,1,2,4,3,1,2,4,3,4,3,2,1,2,3,1,4,3,4,2,1,2,4,1,3,1,3,2,4,2,3,4,1,2,1,3,4,1,4,3,2,3,4,1,2,1,4,2,3,1,3,2,4];
                    toggle = 4;
                case 4
                    temp_sequence = [1,3,4,2,4,2,1,3,4,1,2,3,1,2,3,4,2,4,3,1,4,1,2,3,2,1,3,4,2,1,4,3,4,2,1,3,1,4,2,3,2,1,4,3,2,4,1,3,1,2,4,3,1,2,3,4,3,1,4,2,3,2,1,4,1,3,2,4,1,3,4,2,1,4,3,2,4,1,3,2,1,4,2,3,1,4,2,3,2,3,4,1,2,4,3,1,2,3,1,4,1,3,2,4,2,1,4,3,4,2,1,3,2,1,4,3,1,2,4,3,4,1,2,3,1,2,3,4,1,4,3,2,3,1,4,2,3,2,1,4,1,2,3,4,2,1,3,4,3,2,4,1,2,3,1,4,1,4,3,2];
                    toggle = 5;
                case 5
                    temp_sequence = [4,2,1,3,4,2,1,3,4,3,1,2,3,1,2,4,1,2,4,3,4,1,3,2,1,4,3,2,1,4,2,3,1,2,4,3,1,4,2,3,1,4,3,2,1,4,3,2,3,4,2,1,2,4,1,3,4,2,1,3,2,1,4,3,1,4,2,3,4,3,1,2,3,2,4,1,4,3,2,1,4,1,3,2,4,2,1,3,2,4,1,3,2,3,1,4,3,1,2,4,2,1,4,3,4,2,1,3,2,4,3,1,3,4,1,2,3,2,4,1,3,1,4,2,1,4,2,3,4,3,1,2,3,1,4,2,1,2,3,4,1,2,4,3,4,1,2,3,2,1,4,3,2,4,3,1,4,3,2,1];
                    toggle = 6;
                case 6
                    temp_sequence = [2,4,3,1,3,1,4,2,3,1,4,2,1,4,3,2,1,3,4,2,4,2,1,3,1,3,2,4,3,1,2,4,3,4,2,1,2,1,3,4,2,1,3,4,3,4,1,2,1,3,2,4,1,4,3,2,3,1,2,4,2,4,3,1,3,1,2,4,1,3,4,2,3,4,1,2,3,4,1,2,3,4,2,1,3,4,1,2,3,1,4,2,1,2,3,4,3,1,4,2,3,2,1,4,2,3,1,4,1,2,3,4,1,3,2,4,3,1,4,2,3,1,2,4,3,1,4,2,4,3,2,1,2,3,1,4,3,1,4,2,4,3,2,1,4,1,3,2,3,4,1,2,3,2,1,4,1,3,2,4];
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
            trial = trial + 1;
        end
    end
end
%TRIAL(trial).target=NaN;



HideCursor;
%map of pressed keys
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
%CHANGE THIS FOR DEBUGGING JOYSTICK
screen_properties.joystick_idx = 1;
screen_properties.joystick_sensitivity = 25;
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

%initialize "home" position (origin in this case)
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
position.target_onset = [0];
position.target_movement_onset = [0];
position.home_movement_onset = [0];
position.home_onset = [0];
position.error_x_dproj_vec= [0];
position.dproj_vec = [0];
position.magerror_vec = [0];




%ititialize cursor properties
cursor_subtends.screen_properties = screen_properties;
cursor_subtends.stim_width_mm = 5;
cursor_subtends.stim_height_mm = 5;

trial_type = [0];

i = 1;
draw = 'target';
next = true;
BETWEEN_STIMULI=false;

%setup catch screen
%%%%%%%%%%%%%%%%%%%%
subtends.stim_width_mm = 20;
subtends.stim_height_mm = 20;
subtends.center_on_position=screen_properties.origin + [0 0 0 .2*screen_properties.height_res_pix];
stationary_rect = rect_subtend_distance_mm(subtends);

scale_width = -.2;
while(1)
    [keyPressed sec] = checkKeyboard;
    if(any(keyPressed==KbName('space')))
        break;
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
%%%%%%%%%%%%%%%%%%%


position.t0 = GetSecs();
SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);

while(i<length(TRIAL)+1)
    %disp(num2str(i))
    target_num = TRIAL(i).target;
    %get polar "theta" for current target
    if TRIAL(i).perturb
        position.target_theta = target(target_num).polar_coords(1) - PERTURBATION_ANGLE;
    else
        position.target_theta = target(target_num).polar_coords(1);
    end
    
    %update position, velocity, acceleration, error in heading direction

    
    if TRIAL(i).perturb
        if BETWEEN_STIMULI
            position = center_position(position);
            SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
        else
            position = update_joystick_data_perturb(position,PERTURBATION_ANGLE);
        end
    else
        if BETWEEN_STIMULI
            position = center_position(position);
            SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
        else
            position = update_joystick_data(position);
        end
    end    
    %get polar coords of current pos
    mag_d = position.polar(size(position.polar,1),2);
    theta_d = position.polar(size(position.polar,1),1);

    %draw cursor centered on the current position
    cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
    cursor_rect = rect_subtend_distance_mm(cursor_subtends);
    %Screen('FillOval',mainWin,white,cursor_rect);
    Screen('DrawTexture',mainWin,cursor_tex,[],cursor_rect);
    
    %toggle between drawing targets and home position
    switch draw
        case 'home'
            
            
%            Screen('FrameOval', win, white ,home.rect,4,4);
            
            if next
                BETWEEN_STIMULI=true;
                position.home_onset(i) = GetSecs - position.t0;
                next = false;
            end
            
%            if mag_d < TARGET_RADIUS_MM
            if GetSecs - position.t0 - position.home_onset(i) > WAIT_TIME_POST_TARGET
%                PsychPortAudio('Stop',sound);
                time = GetSecs - position.t0;
                draw = 'target';
                position.time_to_home_from_home_onset(i) = time - position.home_onset(i);
%                position.time_to_home_from_home_movement_onset(i) = time - position.home_movement_onset(i);
                next = true;
                BETWEEN_STIMULI=false;
                
            end
            
        case 'target'
%            Screen('FillOval', mainWin, black ,target(target_num).rect);
            Screen('DrawTexture',mainWin,target_tex,[],target(target_num).rect);
            if next
                position.target_onset(i) = GetSecs - position.t0;
                next = false;
            end
            
            if mag_d > TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM
                
                
                draw = 'home';
                
                if strcmp(TRIAL(i).type,'RD')
                    trial_type = 1;
                elseif strcmp(TRIAL(i).type,'SD')
                    trial_type = 2;
                elseif strcmp(TRIAL(i).type,'RI')
                    trial_type = 3;
                end
                
                time = GetSecs - position.t0; 
                position.home_onset(i) = time;
                position.time_to_target_from_target_onset(i) = time - position.target_onset(i);
                %position.time_to_target_from_target_movement_onset(i) = time - position.target_movement_onset(i);
                
                %update end point error vector
                position.end_point_error(i) = position.error_theta(length(position.error_theta));
                if abs(position.end_point_error(i)) < MAX_SUCCESS_THETA
                   PsychPortAudio('FillBuffer',audiohandle,pleasant_data);
                   PsychPortAudio('Start',audiohandle,1,0,1);
                   feedback = 1;
                else
                   PsychPortAudio('FillBuffer',audiohandle,unpleasant_data);
                   PsychPortAudio('Start',audiohandle,1,0,1);
                   feedback = 0;
                end
                
                
                RESAMPLE_TO = .1;
               
                [interpolated_error,reconstructed_proj_mag] = interpz(RESAMPLE_TO,position.magerror_vec,position.dproj_vec);

                
                
                
                position.max_error_mag(i) = nanmax(interpolated_error);
                max_index = find(interpolated_error==position.max_error_mag(i));
                position.max_error_pos_as_percent(i) = 100*reconstructed_proj_mag(max_index(1))/reconstructed_proj_mag(end);
                position.error_accumulated_over_trajectory(i) = nansum(interpolated_error*RESAMPLE_TO);
                position.avg_error_mag(i) = nanmean(interpolated_error);

                %disp([test2 position.avg_error_mag(i)])
                
                %write data to file to avoid shit getting too damn big
                formatstring = '%d\t%d\t%d\t%d\t%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n';
                fprintf(fileid,formatstring,[repmat(i,length(position.xyt(2:end,3)),1),repmat(trial_type,length(position.xyt(2:end,3)),1),(position.xyt(2:end,3)>position.target_onset(i)),repmat(target_num,length(position.xyt(2:end,3)),1),position.xyt(2:end,3),position.xyt(2:end,1:2), position.velxyt(2:end,1:2), position.accelxyt(2:end,1:2), position.error_vec(2:end,1:2), position.proj_onto_targ(2:end,1:2)]');
                
                
                summaryformatstring = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n';
                fprintf(summaryfile,summaryformatstring,[i,trial_type,target_num,position.time_to_target_from_target_onset(i),position.end_point_error(i),position.avg_error_mag(i),position.error_accumulated_over_trajectory(i),position.max_error_mag(i),position.max_error_pos_as_percent(i),feedback]);
                
                
                %saved
                position.xyt = [0 0 0];
                position.velxyt = [0 0 0];
                position.accelxyt = [0 0 0];
                position.error_vec = [0 0];
                position.proj_onto_targ = [0 0];
                trial_type = [0];
                

                position.error_theta = [0];
                
                %discarded
                position.rect = [0 0 0 0];
                position.polar = [0 0];
                position.error_x_dproj_vec = [0];
                position.dproj_vec = [0];
                position.magerror_vec = [0];

                

                next = true;
                %disp(TRIAL(i).type)
                %disp(i)
                i = i+1;

                
                SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
            end
    end
    

    Screen('Flip',mainWin);
    

    
    
end

%%Handle post test 

stop = KbName('return');
string = '';
disp_char = '';
while(1)
    
    [keyPressed sec]= checkKeyboard;
    if any(keyPressed)
        if (any(keyPressed(1)==stop))
            break;
        elseif any(keyPressed(1) == KbName('space'))
            disp_char = [' '];
        elseif ispc && any(keyPressed(1) == KbName('backspace'))
            string(end) = '';
        elseif ismac && any(keyPressed(1) == KbName('delete'))
            string(end) = '';
        else
            disp_char = KbName(keyPressed(1));
        end
        string = horzcat(string,disp_char);
        disp_char = '';
        keyPressed = 0;
    end
    
    Screen('TextSize',mainWin,25);
    
    DrawFormattedText(mainWin,string,'center','center',black);
  
    DrawFormattedText(mainWin,'Did you notice anything interesting about the golf targets?\nPress "Enter" when finished typing.','center',.2*screen_properties.height_res_pix,black);
    
    Screen('Flip',mainWin);
    
end
postformatstring = '%s\t%s\n';

type_of_post='Guess about targets:';
fprintf(postfile,postformatstring,type_of_post,string);




done_with_post_test = false;
number = [KbName('1') KbName('2') KbName('3') KbName('4') KbName('1!') KbName('2@') KbName('3#') KbName('4$')];
counter = 0;
disp_number = '';
while(counter < length(SEQUENCE))
    
    if done_with_post_test
        break;
    end
    
    [keyPressed sec]= checkKeyboard;
    if (any(keyPressed(1)==number))
        counter = counter + 1;
        disp_char = KbName(keyPressed(1));
        disp_number(counter) = disp_char(1);
        keyPressed = 0;
    end
    
    
    for i = 1:NUMBER_OF_TARGETS
        Screen('TextSize',mainWin,50);
        Screen('DrawTexture',mainWin,target_tex,[],target(i).rect);
        DrawFormattedText(mainWin,num2str(i),target(i).rect(1)+.25*(target(i).rect(3)-target(i).rect(1)),target(i).rect(2),[255 0 0]);
    end
    
    Screen('TextSize',mainWin,25);
    
    DrawFormattedText(mainWin,disp_number,'center','center',black);
  
    DrawFormattedText(mainWin,'During some trials\n you were presented with a \nsequence of locations. \n\nEnter the 8 item sequence:','center',.2*screen_properties.height_res_pix,black);
    

    
    Screen('Flip',mainWin);
    
end
postformatstring = '%s\t%s\n';

type_of_post='numbers:';
disp(disp_number)
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


targ_thetas = [];
for i = 1:length(target)
    targ_thetas(i) = target(i).polar_coords(1);
end
counter = 0;
position.t0_post = GetSecs;
next=true;
draw='home';
disp_number = '';
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
    
    position = update_joystick_data(position);
    
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
    
    position = update_joystick_data(position);
    
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
    
end
disp(disp_number)
type_of_post='forced crct:';
fprintf(postfile,postformatstring,type_of_post,disp_number);





sca
%save('./util/cb.mat','counterbal');




end
