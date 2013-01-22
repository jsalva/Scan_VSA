function testing()
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
%glRotatef(-90, 1,0,0);
%glRotatef(18,0,1,0);
i = 1;
rotheta_last = 0;
rotaxis_x_last = 0;
rotaxis_y_last = 0;
while ~KbCheck
  
  glClear;
  transx = -.01;
  transy = -.01;
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
  Screen('Flip', win);
  Screen('BeginOpenGL', win);

end;
sca
Screen('CloseAll')
end

