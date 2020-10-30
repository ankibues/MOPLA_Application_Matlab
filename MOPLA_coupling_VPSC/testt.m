
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
%{
WK= [1];%,.50,.75,.90,0];     %kk1

A1=[5;5;1]; A2=[5;1;1]; A3=[1;1;1];
AA= [5;5;1];%[A1,A2,A3];          %kk3
 
ANG= generateANG();
 
ANGG= ANG;%ang2,ang3];    %kk4
strain = [7.2,7.4,7.6,7.8];     %kk5 


%*************************************************************************
n=1;
%  generating 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 

 epsilonII=0.5;                                     % initial strain rate invariant  

   for kk1=1:numel(strain)
    

   
   Multi_Parameters_Iso_HEM(1,10,[5;5;1],[330;0;0],strain(kk1),1,1,Jd,b)
         
         
   end
   
   %}
%{
L= [0 1 0;0 0 0;0 0 0];
tincr= .05;
for k=1:1000
    [SI,Gamma,~]=Calc_Gamma(L,tincr,k);
    if SI>=6
        break
    end
end
%}

% plotting vorticity inside the inclusion
figure
plotVorticity(1,0.5,[5;3;1],[45;45;135],6)
hold on
%plotVorticity(1,5,[5;3;1],[45;0;0],6)
hold on
%plotVorticity(1,5,[5;3;1],[135;0;0],6)
hold on
%plotVorticity(1,8,[5;4;1],[0;0;0],6)
hold on
plotVorticity(1,5,[5;4;1],[45;0;0],6)
hold on
plotVorticity(1,2,[5;3;1],[0;90;35],6)
hold on
plotVorticity(.5,2,[5;3;1],[0;45;0],6)
hold on
plotVorticity(1,8,[5;2;1],[0;0;45],6)
%{
figure
plotVorticity(1,8,[5;3;1],[0;90;0],6.8)
hold on
plotVorticity(0,5,[5;1;1],[0;0;60],6.5)
hold on
plotVorticity(.5,5,[5;1;1],[0;0;150],6.5)
hold on
plotVorticity(.8,5,[5;5;1],[0;0;0],6.5)
hold on
plotVorticity(.9,5,[5;5;1],[0;0;150],6.5)

figure
plotVorticityAniso(1,0.5,[1;1;1],[45;45;0],6.5,2)
hold on
plotVorticityAniso(1,2,[5;5;1],[180;90;0],6.5,2)
hold on
plotVorticityAniso(1,5,[5;5;1],[45;90;0],6.5,2)
hold on
plotVorticityAniso(1,8,[5;1;1],[90;90;0],6.5,2)
hold on
plotVorticityAniso(1,.5,[1;1;1],[45;45;0],6.5,5)
%}