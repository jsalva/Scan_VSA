function coords = polar2cartesian(screen_properties,theta,magnitude_mm)

%convert magnitude in mm to x and y components
[x_mm,y_mm] = pol2cart(theta,magnitude_mm);

%convert mm to pix independently for x and y
x_pix = x_mm*screen_properties.width_res_pix/screen_properties.width_mm;
y_pix = y_mm*screen_properties.height_res_pix/screen_properties.height_mm;

coords = [x_pix,y_pix];

return
end