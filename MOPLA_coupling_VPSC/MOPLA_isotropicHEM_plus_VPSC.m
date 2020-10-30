%This script first uses the MOPLA, to determine the partitioned flow in the
%RDE (velocity gradient Tensor L) which is then used as boundary condition
%to the VPSC code.
% Note: this file generates a variable L for number of steps the
% deformation is carried out. Then, this L is put in the VPSC folder in the specified input format,
% and VPSC code is run. Once we get the texture file, it is plotted using
% MTEX toolbox.

%****************************IMPORTANT NOTE********************************

% Need to launch MTXEX toolbox before you run the code
% startup_mtex

%**************************************************************************
%
% Motion of single ellipsoid in isotropic HEM 
% Shear sense investigation and free slip case
addpath(genpath('D:\Ankit\Dropbox\MOPLA_Ankit_Matlab\'));

%  Input parameters--------------------------------------------------------
  Wk= 1;
  gamma=1;
  epsilon= .5*(((1/(Wk^2))-1)^(.5));
 if Wk==0
     gamma=0;
     epsilon=1;
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 else
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 end 
 D = 0.5 * (L + L');
% normalizing with respect to strain rate
 L= L/norm(D);

%  Matrix flow field
   %L     = [0 1 0;0 0 0;0 0 0]; 
   
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 500;
  
%  Power law coefficients for matrix and clast
   Nm=3;
   Nc=1;
   r = .5;
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
   
% for stress invariat for normalization :
   Sigma = 2* Multiply(Jd,D);                                         % far field stress value( sigma= 2*Jd:E)  
   Sigma_Inva = inva(Sigma);

   
% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                          % initial length of the ellipsoid axis
   a= [10;10;1];
   AA=a;
   kkk= 30;
   ang = [kkk;kkk; kkk]  ;                                                     % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                         % initial strain rate invariant
% calculations to follow:
 dotprodcheck= zeros(steps,1);
 r_new=zeros(1,steps);
 L_steps= zeros(3,3,steps);
 D_steps= zeros(3,3,steps);
 WE_steps1=zeros(3,3,steps);
 theta_steps= zeros(3,3,steps);
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
  
for i=1:steps
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q';                                      % w is the vorticity of material, in ellipsoid's frame of reference 
   % [S,~,PI] = SnpIn(a,Jd,Js,Ja);
   
    %****************************
    % Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el,PI_el] = SnP(a);

% Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green
% tensor for pressure)
% p
  p   = zeros(3,3);
  for j=1:3
      p(j,j) = -1/3* (S_el(j,j,1,1)+ S_el(j,j,2,2)+ S_el(j,j,3,3)); 
  end
% S  
  S       = S_el;
  for k=1:3
      for l=1:3
          S(k,k,l,l) = p(k,k)+ S_el(k,k,l,l);
      end
  end
% PI  
  PI      = PI_el;
   %*******************************
    %}
    invS = Inverse(S,b);                    
    E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
    dE= E(1:3,1:3);  % dE is strain rate tensor in ellipsoid
    r_new(1,i)= E(1,4);
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    theta_steps(:,:,i)= theta;
    qanew= QL(a,q,theta,dE,tincr); 
    %updating L_ellipsoid
    ll = dE + we;
    ll= q'*ll*q;
    L_steps(:,:,i)= ll;
    wEE= q'*we*q;
    D_steps(:,:,i)= q'*dE*q;
    WE_steps1(:,:,i)= wEE;
    wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
    dotprod= dot(wvec,zaxis);
    dotprodcheck(i,1)= dotprod; % checking the shear sense !
    
    
    %{
    % giving free slip after 50 steps of deformation
    if i > 50
        
        invSS= Inverse(Jd-S,b);
        we = w + Multiply(Contract(PI,invSS),d);
        wEE= q'*we*q;
       wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
       dotprod= dot(wvec,zaxis);
       dotprodcheck(i,1)= dotprod; % checking the shear sense !
        theta= freeslip(a,w,d,S) ;    % Condition for free slip
        theta_steps(:,:,i)= theta;
        dE= Multiply(invSS,d);
        qanew= QL(a,q,theta,dE,tincr);
        
        ll = dE + we;
        L_steps(:,:,i)= q'*ll*q;
        wEE= q'*we*q;
        wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
        dotprod= dot(wvec,zaxis);
        dotprodcheck(i,1)= dotprod;
    end
      %}    
    %  write updated Q to Q_evl        
       Q_evl(:,:,i)= qanew(1:3,1:3);
    %  write updated a to A_evl
       A_evl(:,i)= qanew(:,4);
    % update new a and q for next step
   %if i < 51
       a= qanew(:,4);
   %end    
    q = qanew(1:3,1:3);
 [SI,Gamma,~]=Calc_Gamma(L,tincr,i);
 if SI>=2
     break
 end
end
steps=i;

numb= ones(steps,1);
tcr= tincr*ones(steps,1);
Lstep= zeros(steps,9);

for k=1:steps
    numb(k)=k;
    kk= L_steps(:,:,k);
    Lstep(k,:)= reshape(kk',1,9);
end

data= [numb Lstep tcr];


%----------------COMING BACK TO THE VPSC FOLDER----------------------------
%---------------------------------------------------------------------------
cd 'D:\Ankit\Dropbox\VPSC_Modeling\Example\pureshear\';
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% code to create .dat file from the L data 
formatspec=  '%-d      %-12.8f      %-12.8f      %-12.8f      %-12.8f      %-12.8f      %-12.8f      %-12.8f      %-12.8f      %-12.8f      %-f\r\n';
formatspec1= '%-d      %-d      %-d      %-d      %-s      %-s      %-s      %-s      %-s      %-s      %-s\r\n';
formatspec2= '%-s      %-s      %-s      %-s      %-s      %-s      %-s      %-s      %-s      %-s      %-f\r\n';

fname = sprintf('L_iso_Wk=%2.1f_r=%2.1f_shape_%d-%d-%d_theta=%d_gamma=%4.2f_SI=%4.2f.dat',Wk,r,AA(1),AA(2),AA(3),kkk,Gamma,SI);
fileid= fopen(fname,'w');
fprintf(fileid,formatspec1,steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' ');
fprintf(fileid,formatspec2,'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', tincr);
for ii=1:steps
       fprintf(fileid,formatspec,data(ii,:));
end
fclose(fileid);


editvpsc7input(fname); %This function edits VPSC7.in file for VPSC input


command= 'D:\Ankit\Dropbox\VPSC_Modeling\Example\pureshear\vpsc7.win64.exe';
status= system(command);

% Plotting the CPO patterns using MTEX toolbox
PlotTexture_AB(fname);

%{
% code to create .xls file from the L data 
firstlin= {steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' '};
seclin= {'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', tincr};
%datanew= {firstlin;seclin;data};
xlswrite('output.xls',firstlin,'sheet1');
xlswrite('output.xls',seclin,'sheet1',sprintf('A%d',2));
xlswrite('output.xls',data,'sheet1',sprintf('A%d',3));
%}
