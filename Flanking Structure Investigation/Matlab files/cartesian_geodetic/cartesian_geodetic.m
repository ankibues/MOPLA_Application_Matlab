function [Geo] = cartesian_geodetic(X, axis )
%
% Converter Cartesian coordinates to Geodetic coordinates 
%on Triaxial Ellipsoid or Biaxial Ellipsoid or Sphere
%
%   (x/a)^2+(y/b)^2+(z/c)^2=1   Triaxial Ellipsoid Equation
%    Cartesian To Geodetic  x y z ==> B L h

% Parameters:
% * X, [x y z]     - Cartesian coordinatesdata, n x 3 matrix or three n x 1 vectors
% * axis,[a; b; c] - ellipsoid radii  [a; b; c],its axes % along [x y z] axes
%  
%                  For Triaxial ellipsoid ,it must be a > b > c
%
%                  For Biaxial ellipsoid ,it must be a = b > c
%
%                  For Sphere ,it must be a = b = c
%
% Output:
% * Geo,[B,L,h]  -  Geodetic coordinates [latitude(deg);longitude(deg);ellipsoidal height(m)]
% 
%
% Author:
% Sebahattin Bektas, 19 Mayis University, Samsun
% sbektas@omu.edu.tr

format long
%ro=180/pi; % converter Degree to radian
eps=0.000001; % three sholder
a=axis(1);b=axis(2);c=axis(3);
   
x=X(1);y=X(2);z=X(3);
%{
if x==0    % this is done to approximate the case when x=0. To avoid getting NaN.
    x=0.00001;
end
%}
%display('a             b                 c                semi axis')
%fprintf(' %15.4f %15.4f %15.4f   \n',a,b,c)


 %display(' C A R T E S  A N     C O O R D  N A T E S         ')   
%display('X                     Y                           Z         ')
% fprintf('  %15.4f %15.4f %15.4f    \n',x,y,z) 
ex2=(a^2-c^2)/a^2; ee2=(a^2-b^2)/a^2;

E=1/a^2;F=1/b^2;G=1/c^2;
xo=a*x/sqrt(x^2+y^2+z^2);
yo=b*y/sqrt(x^2+y^2+z^2);
zo=c*z/sqrt(x^2+y^2+z^2);

for i=1:50
j11=F*yo-(yo-y)*E;
j12=(xo-x)*F-E*xo;

j21=G*zo-(zo-z)*E;
j23=(xo-x)*G-E*xo;

A=[ j11   j12   0 
    j21   0   j23
    2*E*xo    2*F*yo  2*G*zo  ];

sa=(xo-x)*F*yo-(yo-y)*E*xo;
sb=(xo-x)*G*zo-(zo-z)*E*xo;
se=E*xo^2+F*yo^2+G*zo^2-1;
Ab=[ sa  sb  se]';
bil=-A\Ab;
xo=xo+bil(1);
yo=yo+bil(2);
zo=zo+bil(3);

%fprintf('%5d %18.4f %18.4f %18.4f %8.3f %8.3f %8.3f  \n',i,xo,yo,zo,bil(1),bil(2),bil(3))
if max(abs(bil))<eps
    break
end
end
 %display(' C A R T E S  A N     C O O R D  N A T E S  on ELLPSODAL SURFACE       ')   
 %display('Xo                     Yo                           Zo         ')
 %fprintf('  %15.4f %15.4f %15.4f    \n',xo,yo,zo) 
%display(' G E O D E T  C     C O O R D  N A T E S         ')
 fi=atan(zo*(1-ee2)/(1-ex2)/sqrt((1-ee2)^2*xo^2+yo^2));

 l=atan(1/(1-ee2)*yo/xo);

h=sign(z-zo)*sign(zo)*sqrt((x-xo)^2+(y-yo)^2+(z-zo)^2);

Geo=[fi l h];
%display(' Latitude(degree)   Longitude (degree) Ellipsoidal height (m)')
%fprintf(' %18.8f %18.8f %18.3f   \n',fi,l,h) 
end