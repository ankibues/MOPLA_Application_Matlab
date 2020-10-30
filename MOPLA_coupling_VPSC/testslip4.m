% testing memory

% Generating variable velocity gradient tensors for different initial parameters
% Simulating motion of single ellipsoid in isotropic HEM 

 addpath(genpath('F:\Dropbox\MOPLA_Ankit_Matlab\'));

% Different initial parameters

% Wk (Kinematic vorticity number) ranging from 0-1. Wk= 0,.5,.75,.90,1
% Viscosity ratio of RDE to HEM's Nn (since HEM is planar anisotropic) r=.5, 1.5, 2,5,10
% Inital orientation of the RDE given by theta 0,30,60,90
% Initial shape of the RDE  10:10:1, 10:5:1, 10:2:2(Oblate, triaxial and Prolate)

%*************************************************************************
% NOTE: these are the following parameters that are changed in the nested
% for-loops to run simulations for various initial conditions

WK= [.75,.90,0];     %kk1
R=[.5,2,5,8];       %.5,%kk2

A1=[5;3;1]; A2=[5;2;1]; A3=[5;4;1]; A4=[5;1;.999]; A5=[5;5;1]; 
AA= [A1,A2,A3,A4,A5];          %kk3
 
ANG= generateANG();
 
ANGG= ANG;%ang2,ang3];    %kk4
strain = [2,4,6,8];     %kk5 
Nm=1;                 %kk6 
Nn= 1;                %kk7

%*************************************************************************
n=1;
%  generating 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 

 epsilonII=0.5;                                     % initial strain rate invariant  

   for kk1=1:numel(WK)
    

   

         
         for kk2=1:numel(R)

               for kk3=1:(numel(AA)/3)
  
                        for kk4= 1:(numel(ANGG)/3)
                                  
                            for kk5=1:numel(strain)
                                tic
                                 
                                %{
                                if kk1==1 && kk1==2
                                    break
                                end
                                if kk1==3 && kk2<5
                                    break
                                end
                                
                                if kk1<3 && kk2<5 && kk3<3 && kk4<8 && kk5<6
                                    break
                                end
                               %}
                               
                               
                                       
                                Wk=WK(kk1);
                                r= R(kk2);
                                a= AA(:,kk3);
                                ang= ANGG(:,kk4);
                                strI= strain(kk5);
                                
                                Multi_Parameters_Iso_HEM(Wk,r,a,ang,strI,Nm,Nn,Jd,b);  
                               
                               %}
                                disp('A Simulation completed')
                                disp(n)
                                n=n+1;
                                toc
                           
                            end
                        end
               end
         end
   end
   
   


