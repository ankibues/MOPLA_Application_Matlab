a=1;
b=1; 
x= -2:.01:2;
y=-2:.01:2;
Nm= 1;
ang = 0;
q= [cos(ang), sin(ang); -sin(ang), cos(ang)];
E= [1,0;0 -1];
e= q*E*q';
[X,Y] = meshgrid(x,y); 
Pdiff= zeros(sqrt(numel(X)));
for i=1:sqrt(numel(X))
    for j=1:sqrt(numel(X))
        if((X(i,j)^2)/(a^2) + (Y(i,j)^2)/(b^2) > 1)
            rho = (X(i,j)^2 + Y(i,j)^2)^(1/2);
            thetaa= atand(Y(i,j)/X(i,j));
            theta=deg2rad(thetaa);         
            r=2;
            sigma= 2*Nm*E;
            sigma_inva= sigma(1,1)^2;
            h = (a^2 - b^2)/ rho^2;
            R= a/b;
            F= (1+ h^2 - 2*h*cos(2*theta))^(1/2);
            Etilde(1,1) =  e(1,1)* 2*R*(1-r)/(R^2 +2*r*R+1);
            Etilde(1,2) = e(1,2)* (R^2 + 1)*(1-r)/(r*(R^2) + 2*R + r);
            A= ((F + 1 - (h*cos(2*theta)))^(1/2))/(2^(1/2))*F;
            if (theta>0 && theta< pi/2) || (theta > pi && theta< 3*pi/2)
                B= ((F - (1 - (h*cos(2*theta))))^(1/2))/(2^(1/2))*F;
            else
                B= -((F - (1 - (h*cos(2*theta))))^(1/2))/(2^(1/2))*F;
            end      
            Pdiff(i,j) = (8*Nm*R*(R+1)/(R-1))*(((1-A)*(Etilde(1,1)/(2*R)))-(B*(Etilde(1,2)/(R^2 +1))));
            Pdiff(i,j) = Pdiff(i,j)/sigma_inva;
        else
            Pdiff(i,j) = -1.2;
           
        end
    end
end
figure
contourf(X,Y,real(Pdiff))

    
   
  
    
    

    
    %[C, h] = contourf(X,Y,real(Pdiff));
   
   
   % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

   % print -deps graph.eps

 

