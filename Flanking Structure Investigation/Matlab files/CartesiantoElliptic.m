function [epsilon,eta]= CartesiantoElliptic(x,y,A)

a= A(1);  % long semi axis
b= A(2);  % short semi axis
c= (a^2 - b^2)^(1/2);

B= x^2 + y^2 - c^2;

p= (-B + (B^2 + 4*c^2*y^2)^(1/2))/(2*c^2);
q= (-B - (B^2 + 4*c^2*y^2)^(1/2))/(2*c^2);

eta0= asin(p^(1/2));

if x>=0 && y>=0
    eta= eta0;
elseif x<0 && y>=0
    eta= pi- eta0;
elseif x<=0 && y<0
    eta= pi + eta0;
elseif x>0 && y<0
    eta= 2*pi- eta0;
end

epsilon= log(1- 2*q + 2*(q^2-q)^(1/2))/2;

end
