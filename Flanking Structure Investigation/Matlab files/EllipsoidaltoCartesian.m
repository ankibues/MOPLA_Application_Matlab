function [X,Y,Z]= EllipsoidaltoCartesian(lat,long,h,a)

ex= (a(1)^2 - a(3)^2)/a(1)^2;
ee= (a(1)^2 - a(2)^2)/a(1)^2;
N= a(1)/sqrt(1-(ex*(sin(lat))^2)-(ee*(cos(lat))^2*(sin(long))^2));

X= (N+h)*cos(lat)*cos(long);

Y=(N*(1-ee)+h)*cos(lat)*sin(long);

Z=(N*(1-ex)+h)*sin(lat);

end