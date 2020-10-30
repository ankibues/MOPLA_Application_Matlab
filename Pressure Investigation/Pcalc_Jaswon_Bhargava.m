a=4;
b=3;   
[X,Y] = meshgrid(-5:.1:5,-5:.1:5); 

 rho = (X.^2 + Y.^2).^(1/2);
 thetaa= atand(Y./X);
 theta=deg2rad(thetaa);
    ang = pi/3;
    q= [cos(ang), sin(ang); -sin(ang), cos(ang)];
    E= [1,0;0 -1];
    e= q*E*q';
    Nm= 20* rand;
    r=2;
    sigma= 2*Nm*E;
    sigma_inva= sigma(1,1)^2;
    h = (a^2 - b^2)./ rho.^2;
    R= a/b;
    F= (1+ h.^2 - 2*h.*cos(2*theta)).^(1/2);
    Etilde(1,1) =  e(1,1)* 2*R*(1-r)/(R^2 +2*r*R+1);
    Etilde(1,2) = e(1,2)* (R^2 + 1)*(1-r)/(r*(R^2) + 2*R + r);
    A= ((F + 1 - (h.*cos(2*theta))).^(1/2))/(2^(1/2)).*F;
    B= zeros(sqrt(numel(A)));
    for i=1:sqrt(numel(A))
        for j=1:sqrt(numel(A))
            if (theta(i)>0 && theta(i)< pi/2) || (theta(i) > pi && theta(i) < 3*pi/2)
                B(i,j)= ((F(i,j) - (1 - (h(i,j)*cos(2*theta(i,j)))))^(1/2))/(2^(1/2))*F(i,j);
            else
                B(i,j) =  -((F(i,j) - (1 - (h(i,j)*cos(2*theta(i,j)))))^(1/2))/(2^(1/2))*F(i,j);
            end      
        end
    end
       
  
    Pdiff =  (8*Nm*R*(R+1)/(R-1))*(((1-A)*(Etilde(1,1)/(2*R)))-(B*(Etilde(1,2)/(R^2 +1))));
    Pdiff = Pdiff/sigma_inva;

    
    [C, h] = contourf(X,Y,real(Pdiff));
   
   
   % set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

   % print -deps graph.eps

 

