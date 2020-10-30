%load('pressurefield_r=0.1_m=10.mat','pressure');
% note that all the if statements are to make sure that the grid is
% symmetric around the axial planes of the coordinate system so that the
% symmetry of the flow can be utilized to get the pressure/stress invariant
% fields

xgv     = -4:.15:4;
if rem(4,.15)~=0
    xgv= xgv-.05;
    xgv(1)=[];
end

ygv     = -3:.15:3;
if rem(3,.15)~=0
    ygv= ygv-.05;
    ygv(1)=[];
end

 zgv     = -2:.15:2;
 if rem(2,.15)~=0
   zgv= zgv+.05;
 end
 
[X,Y,Z] = meshgrid(xgv,ygv,zgv);
stress = ThreeDfield(stress);
%pressure1=smoothdata(pressure,'sgolay');
figure
%ellipsoidplot(3,1.5,1)
isosurface(X,Y,Z,stress,.2)
isosurface(X,Y,Z,stress,.4)
isosurface(X,Y,Z,stress,1.4)

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