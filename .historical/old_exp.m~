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
            position = update_pos_data_perturb(position,PERTURBATION_ANGLE);
        end
    else
        if BETWEEN_STIMULI
            position = center_position(position);
            SetMouse(screen_properties.origin(3),screen_properties.origin(4), mainWin);
        else
            position = update_pos_data(position);
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
