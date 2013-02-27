function [x,y] = joystick_emulate_mouse(device, sensitivity)
global last_refresh;
%disp(last_refresh)
x_drive = Gamepad('GetAxis',device,1);
if x_drive >=0
    dx = x_drive/abs(AxisRight);
else
    dx = x_drive/abs(AxisLeft);
end

y_drive = Gamepad('GetAxis',device,2);
if y_drive >=0
    dy = y_drive/abs(AxisDown);
else
    dy = y_drive/abs(AxisUp);
end

[x_old,y_old] = GetMouse();
%disp([dx,dy]);

if GetSecs - last_refresh > .01
    x = round(x_old + sensitivity*dx);
    y = round(y_old + sensitivity*dy);
    last_refresh = GetSecs;
else
    x = round(x_old);
    y = round(y_old);
end

SetMouse(x,y);
end