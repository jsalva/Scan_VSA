function position = update_absolute_pos_data(position,reverse_y)
%obtain pixel to distance conversion from position object properties
%(for readability)
global TARGET_DIST_FROM_CENTER_MM
global JOYSTICK_MAGNITUDE_SCALING

x_pix2mm = position.screen_properties.width_mm/position.screen_properties.width_res_pix;
y_pix2mm = position.screen_properties.height_mm/position.screen_properties.height_res_pix;

%origin coords
origin = position.screen_properties.origin;
tmp = Gamepad('GetAxisRawMapping',1,1);
jstick = tmp(1);
x_element = tmp(2);

tmp = Gamepad('GetAxisRawMapping',1,2);
jstick = tmp(1);
y_element = tmp(2);

elements = PsychHID('Elements',jstick);

x_min = elements(x_element).rangeMin;
x_max = elements(x_element).rangeMax;
x_range = x_max - x_min;
x_mid = round(x_range/2);

y_min = elements(y_element).rangeMin;
y_max = elements(y_element).rangeMax;
y_range = y_max - y_min;
y_mid = round(y_range/2);

x_state = PsychHID('RawState',jstick,x_element);
y_state = PsychHID('RawState',jstick,y_element);

x_state_from_center = x_state - x_mid;
y_state_from_center = y_state - y_mid;


[th,mg] = cart2pol(x_state_from_center,y_state_from_center);

th_restricted = mod(th,pi/2);
th_crit = atan2(y_range,x_range);
if th_restricted > th_crit
    y_max = y_max - y_mid;
    x_max = y_max/tan(th_restricted);
elseif th_restricted <= th_crit
    x_max = x_max - x_mid;
    y_max = x_max*tan(th_restricted);
end
[th_max,mg_max] = cart2pol(x_max,y_max);

mg_ratio = (mg/mg_max)*JOYSTICK_MAGNITUDE_SCALING;
current_radius_mm = TARGET_DIST_FROM_CENTER_MM*mg_ratio;
[x_mm,y_mm] = pol2cart(th,current_radius_mm);
x_pix = x_mm/x_pix2mm;
y_pix = y_mm/y_pix2mm;
if reverse_y
   y_pix = -y_pix; 
end

%get current position coords
x_pos_pix = x_pix + origin(3);
y_pos_pix = y_pix + origin(4);


%convert to a psychtoolbox "rect"
position.rect(size(position.rect,1)+1,:) = [0 0 x_pos_pix y_pos_pix];

%convert pixel coordinates to x and y components as vector from orgin
x_pos_mm = (x_pos_pix - origin(3))*x_pix2mm;
y_pos_mm = (origin(4) - y_pos_pix)*y_pix2mm;

%save local copy of current time since start
t = GetSecs()-position.t0; 

%save [x,y,t] for current moment
position.xyt(size(position.xyt,1)+1,:) = [x_pos_mm,y_pos_mm,t];

%differentiate x and y with respect to time (velocity)
dxyt = position.xyt(size(position.xyt,1),:) - position.xyt(size(position.xyt,1)-1,:);
dxdt = dxyt(1)/dxyt(3);
dydt = dxyt(2)/dxyt(3);

%save [x',y',t] for current moment
position.velxyt(size(position.velxyt,1)+1,:) = [dxdt dydt t];

%differentiate x' and y' with respect to time (acceleration)
ddxyt = position.velxyt(size(position.velxyt,1),:) - position.velxyt(size(position.velxyt,1)-1,:);
ddxdt = ddxyt(1)/dxyt(3);
ddydt = ddxyt(2)/dxyt(3);

%save [x'',y'',t] for current moment
position.accelxyt(size(position.accelxyt,1)+1,:) = [ddxdt ddydt t];

%convert position vector from cartesian (x,y) to polar (theta,magnitude)
[theta,mag] = cart2pol(position.xyt(size(position.xyt,1),1),position.xyt(size(position.xyt,1),2));
if theta <= 0
    theta = theta+2*pi();
end
position.polar(size(position.polar,1)+1,:)= [theta,mag];

%save local copy of current position
x_y_coords = [x_pos_mm,y_pos_mm];

%calculate unit vector for current target heading
[unitx,unity] = pol2cart(position.target_theta,1);

%project current position onto target unit vector 
projection = (x_y_coords*[unitx unity]')*[unitx unity];

[null,prev_magproj] = cart2pol(position.proj_onto_targ(size(position.proj_onto_targ,1),1),position.proj_onto_targ(size(position.proj_onto_targ,1),2));
position.proj_onto_targ(size(position.proj_onto_targ,1)+1,:) = projection;

%calculate error vector
error =  x_y_coords - projection;
position.error_vec(size(position.error_vec,1)+1,:) = error;

%take cross product to get sign of error (y)
targ_x_error = cross([unitx unity 0],[position.error_vec(size(position.error_vec,1),:) 0]);
signerror = sign(targ_x_error(3));

%take dot product to get sign of projection(x)
signproj = sign([unitx,unity]*position.proj_onto_targ(size(position.proj_onto_targ,1),:)');

%get magnitude of projection (to get x in "correct"
%direction)
[null,magproj] = cart2pol(projection(1),projection(2));

%get magnitude of error (to get y tangent to "correct" 
%direction
[null,magerror] = cart2pol(error(1),error(2));


%calculate error in theta between current position and target position
%(+ = clockwise error, - = counterclockwise error)
%used to get "endpoint error"
position.error_theta(length(position.error_theta)+1) = atan2((magerror*signerror),(magproj*signproj));
position.dproj_vec(length(position.dproj_vec)+1) = abs(magproj - prev_magproj);
position.magerror_vec(length(position.magerror_vec)+1) = magerror;
position.magproj_vec(length(position.magproj_vec)+1) = magproj;
end
