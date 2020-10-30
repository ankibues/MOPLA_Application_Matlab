%load('pressurefield_r=0.1_m=10.mat','pressure');
xgv     = -4:.1:4; 
ygv     = -3:.1:3;
 zgv     = -2:.1:2;
[X,Y,Z] = meshgrid(xgv,ygv,zgv);
pressure = ThreeDfield(pressure);
pressure=smoothdata(pressure,'sgolay');
figure
ellipsoidplot(3,1.5,1)
isosurface(X,Y,Z,pressure,0)
isosurface(X,Y,Z,pressure,.4)
isosurface(X,Y,Z,pressure,-.4)
%{
load('pressurefield_pureshear_r=0.1_m=10.mat','pressure');
pressure1 = ThreeDfield(pressure);
%pressure1=smoothdata(pressure,'sgolay');
figure
isosurface(X,Y,Z,pressure1,0)
isosurface(X,Y,Z,pressure1,1)
isosurface(X,Y,Z,pressure1,-1)
%}
%pressure1=smoothdata(pressure,'sgolay');