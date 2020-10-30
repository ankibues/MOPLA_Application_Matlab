% Motion of single ellipsoid in Power Law Isotropic Material, with Pressure
% calcutions
  %------------------------------This code snippet is to simply add
  %additional folders and sub-folders to the directory so that the code can be run from Linux terminal 
  k=pwd;
  K= extractBefore(k,"/Pressure Investigation");
  addpath(genpath(K));
  %----------------------------------------------------------------------
%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 1 0;0 0 0;0 0 0]; 
   
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.01;                         % step length for computation
%  total steps of the computation
   steps = 10;
%  number of computation steps between output sets 
   mm=10; 
   
%  Power law coefficients for matrix and clast
   Nm=1;
   Nc=1;
   r= 100;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L'); 
   
   Wvec= [2*W(3,2); 2*W(1,3);2*W(2,1)];  % vorticity vector
   zaxis= [0;0;1];
   dotprod1= dot(Wvec,zaxis);
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 
   
% for stress invariant for normalization :
   Sigma = 2* Multiply(Jd,D);                                         % far field stress value( sigma= 2*Jd:E)  
   Sigma_Inva = inva(Sigma);


% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                          % initial length of the ellipsoid axis
   a= [10;5;1];
   ang = [0; 0; 0]  ;                                                     % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                         % initial strain rate invariant
% calculations to follow:
 Pdev_steps= zeros(steps,1);
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
for i=1:steps
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q';                                      % w is the vorticity of material, in ellipsoid's frame of reference 
    [S,p,PI] = SnpIn(a,Jd,Js,Ja);
    invS= Inverse(S,b);
    h= Contract(2*Jd/Nm, invS);                        
    p_in = R_Multiply(p,h);
    E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    qanew= QL(a,q,theta,dE,tincr);       
    
    %  write updated Q to Q_evl        
        Q_evl(:,:,i)=qanew(1:3,1:3);
    %  write updated a to A_evl
        A_evl(:,i)= qanew(:,4);
        
    Pdev_in= contract1(p_in,u2);
    Pdev_steps(i,1)= Pdev_in/Sigma_Inva;     % normalized with respect to far field stress invariant
    a=qanew(:,4);
    q=qanew(1:3,1:3);
end 
%     compute two spherical angles for three axes
      [a1_ang, a2_ang, a3_ang] = ConvertQ2Angs(Q_evl);
xx=1:steps;

plot(xx,Pdev_steps)
saveas(figure,"test.fig")
set(gcf,'Visible','off');