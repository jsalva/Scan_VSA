function position = update_pos_data_perturb(position,PERTURBATION_ANGLE)
%save local copy of current time since start
t = GetSecs()-position.t0; 

%obtain pixel to distance conversion from position object properties
%(for readability)
x_pix2mm = position.screen_properties.width_mm/position.screen_properties.width_res_pix;
y_pix2mm = position.screen_properties.height_mm/position.screen_properties.height_res_pix;

%origin coords
origin = position.screen_properties.origin;

%get current position coords
[x_pos_pix,y_pos_pix] = GetMouse();

%convert position vector from cartesian (x,y) to polar (theta,magnitude)



%convert to a psychtoolbox "rect"


%convert pixel coordinates to x and y components as vector from orgin
x_pos_mm = (x_pos_pix - origin(3))*x_pix2mm;
y_pos_mm = (origin(4) - y_pos_pix)*y_pix2mm;

[theta,mag] = cart2pol(x_pos_mm,y_pos_mm);

if theta <= 0
    theta = theta+2*pi();
end

theta = theta + PERTURBATION_ANGLE;

position.polar(size(position.polar,1)+1,:)= [theta,mag];

[perturb_x_mm,perturb_y_mm] = pol2cart(theta,mag);
perturb_x_pix = perturb_x_mm/x_pix2mm+origin(3);
perturb_y_pix = origin(4) - perturb_y_mm/y_pix2mm;
        
position.rect(size(position.rect,1)+1,:) = [0 0 perturb_x_pix perturb_y_pix];


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
