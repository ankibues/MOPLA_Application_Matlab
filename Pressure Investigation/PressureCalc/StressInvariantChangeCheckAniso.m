% this calculation shows stress invariant values changing for different
% values of orientation angle of ellipse , ang
L     = [0 2 ;0 0]; 
D = 0.5 * (L + L');
ang = pi/4 ;      
SICheck = zeros(1,10);
Cm    = zeros(2,2,2,2);
 for i=1:10  
   
   % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
   % this calculates Cm for every step, with different ang value    
  q = [cos(ang), sin(ang); -sin(ang), cos(ang)];                      %calculates q for different ang
  
Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + (sin(2*ang))^2/m); Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + (cos(2*ang))^2/m);
Cm(2,2,1,2)= Cm(1,2,2,2);
Cm(2,2,1,1)= Cm(1,1,2,2);
Cm(1,2,1,1)=Cm(1,1,1,2);
Cm(2,1,1,1)= Cm(1,2,1,1);
Cm(1,1,2,1)= Cm(1,1,1,2);
Cm(2,1,2,2)= Cm(1,2,2,2);
Cm(2,2,2,1)= Cm(2,2,1,2);
Cm(2,1,1,2)= Cm(1,2,1,2);
Cm(2,1,2,1)= Cm(2,1,1,2);
Cm(1,2,2,1)= Cm(2,1,2,1);
  
% describe D in the clast's coordinate system 
  D_bar  = q * D * q';
  S_bar  = Multiply2D(Cm, D_bar);                % calculates far field stress in the ellipses coordinate system
  SI     = inva2D(S_bar);                        % stress invariant for far-field for normalization               
  SICheck(1,i)= SI;            % this records SI value for different ang. values
  ang= ang + .1 ;              % increasing ang by .1 rad to see the change 
end