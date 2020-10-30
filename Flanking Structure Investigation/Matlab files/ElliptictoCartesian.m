function [x,y]= ElliptictoCartesian(epsilon,eta,A)

a= A(1);  % long semi axis
b= A(2);  % short semi axis
c= (a^2 - b^2)^(1/2);

x= c*cosh(epsilon)*cos(eta);
y= c*sinh(epsilon)*sin(eta);
end