%function [lat,long,h]= CartesiantoEllipsoidal(XX,YY,ZZ,a)
% transforming from Cartesian coordinate system with coordinates(X,Y,Z) to
% a trixaxial ellipsodial coordinate system with coordinates(lat,long,h)
a=[4;3;2];
XX= 0;
YY=0;
ZZ=4;
% Note that we only calculate for first quadrant(XX>=0;YY>=0, ZZ>=0) and use
% the symmetry of the ellipsoid for calculating points in other quadrants

X= abs(XX);
Y= abs(YY);
Z= abs(ZZ);

% obtaining required constants
ax=a(1);
ay=a(2);
b= a(3);
kx=ax^2-b^2/ax;
ky=ay^2-b^2/ay;
Xbar=X/ax;
Ybar=Y/ay;
Zbar=Z/b;
cx=ax^2/b^2;
cy= ay^2/b^2;
cz= ax^2/ay^2;

% computation of foot point (x,y,z)

if Z==0
    if X==0 && Y==0
        x=0;
        y=0;
        z=b;
    elseif X==0
        x=0;
        y=ay;
        z=0;
    elseif Y==0
        x=ax;
        y=0;
        z=0;
    elseif (ky*X)^2 + (kx*Y)^2 < (kx*ky)^2
        x=ax*X/kx;
        y=ay*Y/ky;
        z= b*sqrt(1-(x/ax)^2-(y/ay)^2);
    else
        f=@(m) (cz*Xbar/(cz+m))^2 + (Ybar/(1+m))^2 -1 ;
        m1= -1+ Ybar;
        m2= -1+ sqrt((cz*Xbar)^2+ Ybar^2);
        M= bisection(f,m1,m2);
        
        x= cz*X/(cz+M);
        y= Y/(1+M);
        z=0;
    end
else
    if X==0 && Y==0
        x=0;
        y=0;
        z=b;
    else
        f=@(m) (cx*Xbar/(cx+m))^2 + (cy*Ybar/(cy+m))^2 + (Zbar/(1+m))^2 -1 ;
        m1= -1+ Zbar;
        m2= -1+ sqrt((cx*Xbar)^2+(cy*Ybar)^2+ Zbar^2);
        M= bisection(f,m1,m2);
        
        x= cx*X/(cx+M);
        y= cy*Y/(cy+M);
        
        if M<0 && (ky*X)^2 + (kx*Y)^2 < (kx*ky)^2
            z= b*sqrt(1-(x/ax)^2-(y/ay)^2);
        else
            z= Z/(1+M);
        end
    end
end

ex2= (ax^2-b^2)/ax^2;
ee2= (ax^2-ay^2)/ax^2;


if (1-ee2)*z <= (1-ex2)*sqrt(y^2 + x^2*(1-ee2)^2)
    lat= atan((1-ee2)*z/(1-ex2)/sqrt(y^2 + x^2*(1-ee2)^2));
else
    lat= pi/2 - atan((1-ex2)*sqrt(y^2 + x^2*(1-ee2)^2)/(1-ee2)/z);
end


if y<= (1-ee2)*x
    long= 2*atan(y/((1-ee2)*x + sqrt(y^2 + x^2*(1-ee2)^2)));
else
    long= pi/2 - 2*atan((1-ee2)*x/(y+ sqrt(y^2 + x^2*(1-ee2)^2)));
end

if x==y && x==0
    long=pi/2;
end
h= sqrt((X-x)^2 + (Y-y)^2 + (Z-z)^2);

if (X+Y+Z)<(x+y+z)
    h=-h;
end

% adjusting according to different quadrants

if ZZ<0
    lat=-lat;
end

if XX<0 && YY>=0
    long=pi-long;
elseif XX<=0 && YY<0
    long=long- pi;
elseif XX>0 && YY<0
    long=-long;
end
   

%