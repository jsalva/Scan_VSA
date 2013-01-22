function Run_VSA(subject_id)

addpath ./util
addpath ./data

filename = ['./data/',subject_id,'_data.txt'];
fileid = fopen(filename,'w');

headerstring = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fileid,headerstring,'Trial','Type','Target/Home(1,0)','Target #','Time(s)','Xpos','Ypos','Xvel','Yvel','Xaccel','Yaccel','Xerror','Yerror','Xprojection','Yprojection');
                
%setup experiment paramters
NUMBER_OF_TARGETS = 4;
PERTURBATION_ANGLE = pi()/4;
LENGTH_OF_SEQUENCE = 8;
TARGET_DIST_FROM_CENTER_MM = 100;
TARGET_RADIUS_MM = 5;
WAIT_TIME_POST_TARGET = .1;
cb = load('./util/cb.mat');
counterbal = cb.counterbal;

%sound stuff
InitializePsychSound;


%OpenGL stuff
AssertOpenGL;
screenid=max(Screen('Screens'));
InitializeMatlabOpenGL([],[],[], 0);
multiSample=0;
imagingPipeline=0;
[win , winRect] = Screen('OpenWindow', screenid, 0,[],[],[],0,multiSample, imagingPipeline);
Screen('BeginOpenGL', win);
ar=winRect(4)/winRect(3);
glColor3f(1,1,0);
glEnable(GL.LIGHTING);
glEnable(GL.LIGHT0);
glEnable(GL.DEPTH_TEST);
glMatrixMode(GL.PROJECTION);
glLoadIdentity;
gluPerspective(25,1/ar,0.1,100);
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;
glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);
gluLookAt(3,3,5,0,0,0,0,1,0);
glClearColor(0,0,0,0);
glClear;
Screen('EndOpenGL', win);
Screen('Flip', win);
myimg = imread([PsychtoolboxRoot 'PsychDemos/OpenGL4MatlabDemos/earth_512by256.jpg']);
mytex = Screen('MakeTexture', win, myimg, [], 1);
[gltex, gltextarget] = Screen('GetOpenGLTexture', win, mytex);
Screen('BeginOpenGL', win);
glEnable(gltextarget);
glBindTexture(gltextarget, gltex);
glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.MODULATE);
glTexParameteri(gltextarget, GL.TEXTURE_WRAP_S, GL.REPEAT);
glTexParameteri(gltextarget, GL.TEXTURE_WRAP_T, GL.REPEAT);
glTexParameteri(gltextarget, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
glTexParameteri(gltextarget, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
glMaterialfv(GL.FRONT_AND_BACK,GL.AMBIENT, [ 1 1 1 1 ]);
glMaterialfv(GL.FRONT_AND_BACK,GL.DIFFUSE, [ 1 1 1 1 ]);
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;
gluLookAt(0,0,5,0,0,0,0,1,0);
mysphere = gluNewQuadric;
gluQuadricTexture(mysphere, GL.TRUE); 
i = 1;
rotheta_last = 0;
rotaxis_x_last = 0;
rotaxis_y_last = 0;
Screen('EndOpenGL', win);

[data,freq,bits] = wavread('./util/golfsound.wav');
channels = size(data,2);
soundhandle=PsychPortAudio('Open',[],[],0,freq,channels,[],0.1);
PsychPortAudio('FillBuffer',soundhandle,data');

if counterbal == 1
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [5 10 5 10 5 20 5];
    BLOCK_STRUCTURE = {'RD','SD','RD','SD','RD','RI','RD'};
    counterbal = 2;
elseif counterbal ==2 
    NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK = [5 20 5 10 5 10 5];
    BLOCK_STRUCTURE = {'RD','RI','RD','SD','RD','SD','RD'};
    counterbal = 1;
end
RAND = false;

while (~exist('SEQUENCE','var')) || SEQUENCE(1) == SEQUENCE(end)
    if RAND
        SEQUENCE = sequence_randomizer(LENGTH_OF_SEQUENCE,NUMBER_OF_TARGETS);
    else
        SEQUENCE = [2 3 1 4 3 2 4 1];% 3 4 2 1];
    end
end
trial = 1;
toggle = 1;
for block = 1:length(BLOCK_STRUCTURE)
    %determine if sequence vs not
    if strcmp(BLOCK_STRUCTURE{block},'SD')
        temp_sequence = repmat(SEQUENCE,1,NUMBER_OF_SEQUENCE_REPEATS_PER_BLOCK(block));
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
                case 1
                    temp_sequence = [4,1,3,2,3,2,4,1,2,1,3,4,2,3,4,1,2,3,4,1,4,3,1,2,4,2,1,3,2,3,4,1,2,4,1,3,4,2,3,1,4,1,3,2,4,3,1,2,1,4,3,2,3,4,2,1,4,2,1,3];
                    toggle = 2;
                case 2
                    temp_sequence = [2,4,3,1,4,2,3,1,3,4,1,2,1,2,4,3,1,2,4,3,1,4,3,2,1,4,3,2,4,1,2,3,2,3,4,1,1,4,2,3,4,3,1,2,2,1,4,3,2,4,1,3,1,2,4,3,1,2,3,4];
                    toggle = 3;
                case 3
                    temp_sequence = [2,1,4,3,2,1,3,4,3,4,2,1,4,3,1,2,3,2,4,1,3,1,2,4,3,4,1,2,1,2,4,3,2,4,1,3,4,2,3,1,4,1,2,3,2,4,1,3,2,3,1,4,1,4,2,3,4,3,2,1];
                    toggle = 4;
                case 4
                    temp_sequence = [4,3,2,1,2,1,3,4,3,2,1,4,3,2,4,1,4,2,1,3,1,4,3,2,4,2,1,3,4,3,2,1,2,3,4,1,4,1,3,2,4,3,2,1,1,2,4,3,4,2,1,3,1,3,2,4,4,3,2,1];
                    toggle = NaN;
            end
        end

        for i = 1:length(temp_sequence)
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
TRIAL(trial).target=NaN;


HideCursor
%map of pressed keys
global downmap
[keyIsDown sec keyCode] = KbCheckM;
downmap = zeros(size(keyCode));

% Open win on display monitor
%screenNum = max(Screen('Screens'));
%win = Screen('OpenWindow', screenNum);
%winRect = Screen('Rect', win);

%initialize colors
black = BlackIndex(win);
white = WhiteIndex(win);

%create main screen
Screen('FillRect', win, black);

%setup stimulus target array
delta_theta = 2*pi()/NUMBER_OF_TARGETS;

%specify screen properties
screen_properties.width_mm = winRect(3)/4;
screen_properties.height_mm = winRect(4)/4;
screen_properties.width_res_pix = winRect(3);
screen_properties.height_res_pix = winRect(4);
screen_properties.subj_distance_mm = 300;
screen_properties.origin = winRect/2;
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
SetMouse(screen_properties.origin(3),screen_properties.origin(4), win);
position.rect = [0 0 0 0];
position.xyt = [0 0 0];
position.polar = [0 0];
position.t0 = GetSecs();
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



%ititialize cursor properties
cursor_subtends.screen_properties = screen_properties;
cursor_subtends.stim_width_mm = 5;
cursor_subtends.stim_height_mm = 5;

trial_type = [0];

i = 1;
draw = 'target';
next = true;
BETWEEN_STIMULI = false;
while(1)
    %get current target
    target_num = TRIAL(i).target;%SEQUENCE(mod(i-1,12)+1);
    if isnan(target_num)
        break;
    end
    
    
    %get polar "theta" for current target
    if TRIAL(i).perturb
        position.target_theta = target(target_num).polar_coords(1) - PERTURBATION_ANGLE;
    else
        position.target_theta = target(target_num).polar_coords(1);
    end
    
    %update position, velocity, acceleration, error in heading direction
    if TRIAL(i).perturb
        position = update_pos_data_perturb(position,PERTURBATION_ANGLE);
    elseif BETWEEN_STIMULI
        %saved
        position.xyt = [0 0 0];
        position.velxyt = [0 0 0];
        position.accelxyt = [0 0 0];
        position.error_vec = [0 0];
        position.proj_onto_targ = [0 0];
        trial_type = [0];
                
        %discarded
        position.rect = [0 0 0 0];
        position.polar = [0 0];
        position.error_theta = [0];
    else
        position = update_pos_data(position);
    end    
    %get polar coords of current pos
    mag_d = position.polar(size(position.polar,1),2);
    theta_d = position.polar(size(position.polar,1),1);

    %draw cursor centered on the current position
    cursor_subtends.center_on_position = position.rect(size(position.rect,1),:);
    cursor_rect = rect_subtend_distance_mm(cursor_subtends);
    Screen('FillOval',win,white,cursor_rect);
    Screen('BeginOpenGL', win);
    glClear;
    vec = position.xyt(size(position.xyt,1),1:2) - position.xyt(size(position.xyt,1)-1,1:2);
    %dtheta = theta_d - position.polar(size(position.polar,1)-1,1);
    %[transx,transy] = pol2cart(dtheta,dmag);
    transx = vec(1);
    transy = vec(2);
    scale=.008;
    transx = scale*transx;
    transy = scale*transy;
    radius=.1;
    arclength=norm([transx,transy]);
    rotheta = (arclength/radius)*360/(2*pi);
    [rotaxis_x,rotaxis_y]=pol2cart(atan2(transy,transx)+pi()/2,1);
    unrotangle = -rotheta_last;
    glRotated(unrotangle,rotaxis_x_last,rotaxis_y_last,0);
    glTranslatef(transx,transy,0)
    rotangle = rotheta_last + rotheta;
    glRotated(rotangle, rotaxis_x, rotaxis_y, 0);
    rotaxis_x_last = rotaxis_x;
    rotaxis_y_last = rotaxis_y;
    rotheta_last = mod(rotheta_last + rotheta,360);
    %  glTranslatef(-point(1),-point(2),-point(3))
    %  glRotated(1,1,0,0)
    gluSphere(mysphere, radius, 100, 100);
    Screen('EndOpenGL', win);
    
 
    
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
                PsychPortAudio('Stop',sound);
                time = GetSecs - position.t0;
                draw = 'target';
                position.time_to_home_from_home_onset(i) = time - position.home_onset(i);
%                position.time_to_home_from_home_movement_onset(i) = time - position.home_movement_onset(i);
                next = true;
                BETWEEN_STIMULI=false;
                
            end
            
        case 'target'
            Screen('FrameOval', win, white ,target(target_num).rect,4,4);
            if next
                position.target_onset(i) = GetSecs - position.t0;
                next = false;
            end
            
            if mag_d > TARGET_DIST_FROM_CENTER_MM - TARGET_RADIUS_MM
                PsychPortAudio('Start',soundhandle,1,0,1);
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
                disp(position.end_point_error(i));
                
                %write data to file to avoid shit getting too damn big
                formatstring = '%d\t%d\t%d\t%d\t%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n';
                fprintf(fileid,formatstring,[repmat(i,length(position.xyt(2:end,3)),1),repmat(trial_type,length(position.xyt(2:end,3)),1),(position.xyt(2:end,3)>position.target_onset(i)),repmat(target_num,length(position.xyt(2:end,3)),1),position.xyt(2:end,3),position.xyt(2:end,1:2), position.velxyt(2:end,1:2), position.accelxyt(2:end,1:2), position.error_vec(2:end,1:2), position.proj_onto_targ(2:end,1:2)]');
                
                %saved
                position.xyt = [0 0 0];
                position.velxyt = [0 0 0];
                position.accelxyt = [0 0 0];
                position.error_vec = [0 0];
                position.proj_onto_targ = [0 0];
                trial_type = [0];
                
                %discarded
                position.rect = [0 0 0 0];
                position.polar = [0 0];
                position.error_theta = [0];
                

                next = true;
                i = i+1;
            end
    end
    

    Screen('Flip',win);
    

    
    
end
sca
save('./util/cb.mat','counterbal');

PsychPortAudio('Close',sound)


end