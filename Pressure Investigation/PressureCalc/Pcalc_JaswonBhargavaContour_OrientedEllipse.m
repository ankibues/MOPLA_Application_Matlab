a=3;
b=1;
R= a/b;
x= -5:.05:5;
y=-5:.05:5;
Nm= 1;
 ang = 0;
 L= [2,0;0,-2];
 q= [cos(ang), sin(ang); -sin(ang), cos(ang)]; % I have written the value of Q (Transformation matrix) in 2D. Q= [a1' a2'] (') is for Transpose.(t;?????
E= .5* (L+L');
 e= q*E*q'; % (') is for Transpose. This is for bringing E to ellipses axis.I don't understand this
[X,Y] = meshgrid(x,y); 
Pdiff= zeros(sqrt(numel(X)));
 r=100;
 sigma= 2*Nm*E;
 sigma_inva= inva2D(sigma);
 Pinside = 2*Nm*e(1,1)*(1-r)*(R^2-1)/(R^2+2*r*R+1)/sigma_inva;
 
for i=1:sqrt(numel(X))
    for j=1:sqrt(numel(X))
        if(((X(i,j)*cos(ang)+ Y(i,j)*sin(ang))/a)^2 + ((X(i,j)*sin(ang)- Y(i,j)*cos(ang))/b)^2 > 1)
            rho = (X(i,j)^2 + Y(i,j)^2)^(1/2);
            thetaa= atand(Y(i,j)/X(i,j));
            theta=deg2rad(thetaa)  - ang ;         
            h = (a^2 - b^2)/ rho^2;
            F= (1+ h^2 - 2*h*cos(2*theta))^(1/2);
            Etilde(1,1) =  E(1,1)* 2*R*(1-r)/(R^2 +2*r*R+1);
            Etilde(1,2) = E(1,2)* (R^2 + 1)*(1-r)/(r*(R^2) + 2*R + r);
            A= ((F + 1 - (h*cos(2*theta)))^(1/2))/(2^(1/2))/F;
            if (theta>0 && theta< pi/2) || (theta > pi && theta< 3*pi/2)
                B= ((F - (1 - (h*cos(2*theta))))^(1/2))/(2^(1/2))/F;
            else
                B= -((F - (1 - (h*cos(2*theta))))^(1/2))/(2^(1/2))/F;
            end      
            Pdiff(i,j) = -(8*Nm*R*(R+1)/(R-1))*(((1-A)*(Etilde(1,1)/(2*R)))-(B*(Etilde(1,2)/(R^2 +1))));
            Pdiff(i,j) = Pdiff(i,j)/sigma_inva;
        else
            Pdiff(i,j) = Pinside;  % this 
           
        end
    end
end
figure
contourf(X,Y,real(Pdiff))

    
   
  
    
 

