

% Generating variable velocity gradient tensors for different initial parameters
% Simulating motion of single ellipsoid in isotropic HEM 

 addpath(genpath('D:\Ankit\Dropbox\MOPLA_Ankit_Matlab\'));

% Different initial parameters

% Wk (Kinematic vorticity number) ranging from 0-1. Wk= 0,.5,.75,.90,1
% Viscosity ratio of RDE to HEM's Nn (since HEM is planar anisotropic) r=.5, 1.5, 2,5,10
% Inital orientation of the RDE given by theta 0,30,60,90
% Initial shape of the RDE  10:10:1, 10:5:1, 10:2:2(Oblate, triaxial and Prolate)

%*************************************************************************
% NOTE: these are the following parameters that are changed in the nested
% for-loops to run simulations for various initial conditions

WK= [1];%,.50,.75,.90,0];     %kk1

R=[.5,2,5,8,10];       %kk2

A1=[5;5;1]; A2=[5;1;1]; A3=[1;1;1]; 
AA= [A1,A2,A3];          %kk3
 
ANG= generateANG2();
 
ANGG= ANG;%ang2,ang3];    %kk4
strain = [6.5,7,7.5];     %kk5
MM= [2,5];         
Nm=1;
Nn=1;
%*************************************************************************
n=1;
 %  generating 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 

   % gaussian points calculation   
   gp                = 20;
   [p, w1]            = Gauss(gp);
   ww                = w1 * w1';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
   
   %}

   for kk1=1:numel(WK)
    

   

         
         for kk2=1:numel(R)

               for kk3=1:(numel(AA)/3)
  
                        for kk4= 1:(numel(ANGG)/3)
                                    
                            for kk5=1:numel(strain)
                                for kk6=1:numel(MM)
                                    
                                    tic
                                    
                                    Wk=WK(kk1);
                                    r= R(kk2);
                                    a= AA(:,kk3);
                                    ang= ANGG(:,kk4);
                                    strI= strain(kk5);
                                    %
                                    m= MM(kk6);
                                    Multi_Parameters_Aniso_HEM(Wk,r,a,ang,strI,m,Jd,Ja,Js,b,Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww);  
                                    disp('Simulation no. completed')
                                    %}
                                    %Multi_Parameters_Iso_HEM2(Wk,r,a,ang,strI,Nm,Nn,Jd,b);  
                                
                                    disp(n)
                                    n=n+1;
                                    
                                    %clearvars Wk r a ang strI m
                                end
                           
                            end
                        end
               end
         end
   end
   
   


