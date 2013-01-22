function rect = rect_subtend_vis_angle(subtends_obj)

angle_x_rad = subtends_obj.stim_width_rad;
angle_y_rad = subtends_obj.stim_height_rad;
screen_width_mm = subtends_obj.screen_properties.width_mm;
screen_height_mm = subtends_obj.screen_properties.height_mm;
center_x = subtends_obj.center_on_position(1);
center_y = subtends_obj.center_on_position(2);
x = subtends_obj.screen_properties.width_res_pix;
y = subtends_obj.screen_properties.height_res_pix;
distance_mm = subtends_obj.screen_properties.subj_distance_mm;

x_dist_pix = tan(angle_x_rad/2)*distance_mm*x/screen_width_mm;
y_dist_pix = tan(angle_y_rad/2)*distance_mm*y/screen_height_mm;

rect(1) = (center_x - x_dist_pix);
rect(2) = (center_y - y_dist_pix);
rect(3) = (center_x + x_dist_pix);
rect(4) = (center_y + y_dist_pix);
disp(rect)

end