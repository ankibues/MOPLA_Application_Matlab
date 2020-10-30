function plotVorticity(Wk,r,a,ANG,strI)
% Motion of single ellipsoid in isotropic HEM 

%addpath(genpath('D:\Ankit\Dropbox\MOPLA_Ankit_Matlab\'));
%  Input parameters--------------------------------------------------------
% Wk=1;
 %r=2;
 %a=[5;2;1];
 %ANG=[0;0;0];
% strI=6;
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
 %L= L/norm(D);

%  Matrix flow field
  % L     = [0 1 0;0 0 0;0 0 0]; 
   
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 300;
  
%  Power law coefficients for matrix and clast
   Nm=1;
   Nc=1;
   
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L'); 
   Wvec= [2*W(3,2); 2*W(1,3);2*W(2,1)];  % macroscale vorticity vector
   mag_Wvec= norm(Wvec);
  % zaxis= [0;0;1];
   %dotprod1= dot(Wvec,zaxis);
%  generate 4th-order identity tensors   
   [Jd, ~, ~,~, b] = Jnb(); 
   
% for stress invariat for normalization :
   %Sigma = 2* Multiply(Jd,D);                                         % far field stress value( sigma= 2*Jd:E)  
 

   
% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                          % initial length of the ellipsoid axis
   
  % AA=a;
   %kkk= 60*pi/180;
   %ang = ANG*pi/180;                                                % initial angles of the ellipsoid axis
   q = Q(ANG);
   epsilonII=0.5;                                                         % initial strain rate invariant
% calculations to follow:
strainI= zeros(steps,1);
 dotprodcheck= zeros(steps,1);
 dotPROD=zeros(steps,1);
 r_new=zeros(1,steps);
 L_steps= zeros(3,3,steps);
 D_steps= zeros(3,3,steps);
 WE_steps1=zeros(3,3,steps);
 theta_steps= zeros(3,3,steps);
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
  mag_wvec= zeros(1,steps);
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
    D_steps(:,:,i)= dE;%q'*dE*q;
    WE_steps1(:,:,i)= wEE;
    wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector of the RDE
    mag_wvec(i)= norm(wvec);
    dotPROD(i)= dot(wvec,Wvec);   % this is the dot product of two vorticity vectors - normalized with respect to magnitude
    
    dotprodcheck(i,1)= dotPROD(i)/(mag_wvec(i)*mag_Wvec); % checking the shear sense !
    
    
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
 [SI,~,~]=Calc_Gamma(L,tincr,i);
 strainI(i)= SI;
 if SI>=strI
     break
 end
end
steps=i;
%
numb= ones(steps,1);
%tcr= tincr*ones(steps,1);
Lstep= zeros(steps,9);

for k=1:steps
    numb(k)=k;
    kk= L_steps(:,:,k);
    Lstep(k,:)= reshape(kk',1,9);
end
%DT= repmat(dotprod1,steps,1);

DOTprod= dotprodcheck(1:steps);
%
DOTprod= DOTprod(1:i);
strainI=strainI(1:i);
if r<1
    plot(strainI,DOTprod,'b')
else
    plot(strainI,DOTprod,'r')
end
hold on
%plot(1:steps,DT,'--')
set(gcf,'Visible','on');
end
%}