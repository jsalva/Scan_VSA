function rect = rect_subtend_distance_mm(subtends_obj)

width_mm = subtends_obj.stim_width_mm;
height_mm = subtends_obj.stim_height_mm;
screen_width_mm = subtends_obj.screen_properties.width_mm;
screen_height_mm = subtends_obj.screen_properties.height_mm;
center_x = subtends_obj.center_on_position(3);
center_y = subtends_obj.center_on_position(4);
x = subtends_obj.screen_properties.width_res_pix;
y = subtends_obj.screen_properties.height_res_pix;
distance_mm = subtends_obj.screen_properties.subj_distance_mm;

x_dist_pix = width_mm*x/screen_width_mm;
y_dist_pix = height_mm*y/screen_height_mm;

rect(1) = (center_x - x_dist_pix);
rect(2) = (center_y - y_dist_pix);
rect(3) = (center_x + x_dist_pix);
rect(4) = (center_y + y_dist_pix);

if rect(1) < 0
    rect(1) = 0;
    rect(3) = 2*x_dist_pix;
end
if rect(2) < 0 
    rect(2) = 0;
    rect(4) = 2*y_dist_pix;
end
if rect(3) > x
    rect(1) = x - 2*x_dist_pix;
    rect(3) = x;
end
if rect(4) > y
    rect(2) = y - 2*y_dist_pix;
    rect(4) = y;
end

end
