% Fletcher (2009) result to show how pressure value can differ for pure and simple shear flow.
ep=0;
v= 1;
%ellipse axial ratio
R= 5;
%viscosity ratio of inclusion to that of matrix(Nn)
r = 100; 
Nn =1;
m=25;
sqrtm= m^.5;
% stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  Cm    = zeros(2,2,2,2);
  pdiff_anal= zeros(180,1);
  for i=0:180  
      ang = (pi/180)*i;
      Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + ((sin(2*ang))^2)/m);
      Cm(2,2,2,2)=  Cm(1,1,1,1);
      Cm(1,1,2,2)= -Cm(1,1,1,1);
      Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
      Cm(1,2,2,2)= -Cm(1,1,1,2);
      Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + ((cos(2*ang))^2)/m);
      Cm(2,2,1,2)= Cm(1,2,2,2);
      Cm(2,2,1,1)= Cm(1,1,2,2);
      Cm(1,2,1,1)= Cm(1,1,1,2);
      Cm(2,1,1,1)= Cm(1,2,1,1);
      Cm(1,1,2,1)= Cm(1,1,1,2);
      Cm(2,1,2,2)= Cm(1,2,2,2);
      Cm(2,2,2,1)= Cm(2,2,1,2);
      Cm(2,1,1,2)= Cm(1,2,1,2);
      Cm(2,1,2,1)= Cm(2,1,1,2);
      Cm(1,2,2,1)= Cm(2,1,2,1);

% analytical results from Fletcher 2009  (as in Jiang 2016)
  
    E11= ep*cos(2*ang) + v*sin(2*ang)/2;
    E11m= ep*cos(2*ang) + (v/m)*sin(2*ang)/2;
    E12= -ep*sin(2*ang) + v*cos(2*ang)/2;
    E12m= -ep*sin(2*ang) + (v/m)*cos(2*ang)/2;
    W12= v/2;
    E11bar= -r*E11 + E11m;
    E12bar= -r*E12 + E12m;
    e11analytic = 2*sqrtm*R*E11bar/ (R^2 + 2*r*sqrtm*R + 1);
    pdiff_analytic = 2*Nn*(R^2 - 1)*(E11bar)/(R^2 + 2*r*sqrtm*R + 1);
    E_anal= [E11, E12; E12, -E11];
   
    S_anal=  Multiply2D(Cm, E_anal);                            
    SI_anal= inva2D(S_anal);                      % Analytical stress invariant for far-field for normalization
   
    pdiff_anal(i+1)= pdiff_analytic/SI_anal;
    
  end
  x=0:180;
plot (x,pdiff_anal,'--')
xlabel('Angle in degrees');
ylabel('alpha1 analytical');
   