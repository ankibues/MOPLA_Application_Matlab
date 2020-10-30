%Trying the variable bulk flow field in case of multiscale structural
%modelling
%************** first bulk flow field
%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 1 0;0 0 0;0 0 0]; 
   
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 200;
  
%  Power law coefficients for matrix and clast
   Nm=3;
   Nc=1;
   r= 2;
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
   a= [10;5;1];
   k= 0;
   ang = [0; 0; k*pi/180]  ;                                                     % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                         % initial strain rate invariant
% calculations to follow:
 dotprodcheck= zeros(steps*2,1);
 r_eff= zeros(1,steps);
 L_steps= zeros(3,3,steps*2);
  Q_evl= zeros(3,3,steps*2);
  A_evl= zeros(3,steps*2);
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
    r_eff(1,i)= E(1,4);
    dE= E(1:3,1:3);  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    qanew= QL(a,q,theta,dE,tincr); 
    %updating L_ellipsoid
    ll = dE + we;
    ll= q'*ll*q;
    L_steps(:,:,i)= ll;
    wEE= q'*we*q;
    wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
    dotprod= dot(wvec,zaxis);
    dotprodcheck(i,1)= dotprod; % checking the shear sense !
    
    
    %{
    % giving free slip after 50 steps of deformation
    if i > 30
        
        invSS= Inverse(Jd-S,b);
        we = w + Multiply(Contract(PI,invSS),d);
        wEE= q'*we*q;
       wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
       dotprod= dot(wvec,zaxis);
       dotprodcheck(i,1)= dotprod; % checking the shear sense !
        theta= freeslip(a,w,d,S) ;    % Condition for free slip
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
   %if i < 31
       a= qanew(:,4);
   %end    
    q = qanew(1:3,1:3);
end

%{
%************** second bulk flow field
%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 sqrt(2) 0;0 0 0;0 sqrt(2) 0]; 
   
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L'); 
   Wvec= [2*W(3,2); 2*W(1,3);2*W(2,1)];  % vorticity vector
   zaxis= [0;0;1];
   dotprod1_= dot(Wvec,zaxis);
   
% for stress invariat for normalization :
   Sigma = 2* Multiply(Jd,D);                                         % far field stress value( sigma= 2*Jd:E)  
   Sigma_Inva = inva(Sigma);

% initial state of the ellipsoid (we get from previous simulation) (a and q)
   
                                                    
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
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    qanew= QL(a,q,theta,dE,tincr); 
    %updating L_ellipsoid
    ll = dE + we;
    ll= q'*ll*q;
    L_steps(:,:,i+steps)= ll;
    wEE= q'*we*q;
    wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
    dotprod= dot(wvec,zaxis);
    dotprodcheck(i,1)= dotprod; % checking the shear sense !
    
    
    %{
    % giving free slip after 50 steps of deformation
    if i > 30
        
        invSS= Inverse(Jd-S,b);
        we = w + Multiply(Contract(PI,invSS),d);
        wEE= q'*we*q;
       wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
       dotprod= dot(wvec,zaxis);
       dotprodcheck(i,1)= dotprod; % checking the shear sense !
        theta= freeslip(a,w,d,S) ;    % Condition for free slip
        dE= Multiply(invSS,d);
        qanew= QL(a,q,theta,dE,tincr);
        
        ll = dE + we;
        L_steps(:,:,i)= q'*ll*q;
        wEE= q'*we*q;
        wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
        dotprod= dot(wvec,zaxis);
        dotprodcheck(i,1)= dotprod;
    end
      }    
    %  write updated Q to Q_evl        
       Q_evl(:,:,i+steps)= qanew(1:3,1:3);
    %  write updated a to A_evl
       A_evl(:,i+steps)= qanew(:,4);
    % update new a and q for next step
   %if i < 31
       a= qanew(:,4);
   %end    
    q = qanew(1:3,1:3);
end
%}
%%savinf to excel file

numb= ones(2*steps,1);
tcr= tincr*ones(2*steps,1);
Lstep= zeros(2*steps,9);

for k=1:2*steps
    numb(k)=k;
    kk= L_steps(:,:,k);
    Lstep(k,:)= reshape(kk',1,9);
end
%{
data= [numb Lstep tcr];

firstlin= {2*steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' '};
seclin= {'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', tincr};
%datanew= {firstlin;seclin;data};
xlswrite('output.xls',firstlin,'sheet1');
xlswrite('output.xls',seclin,'sheet1',sprintf('A%d',2));
xlswrite('output.xls',data,'sheet1',sprintf('A%d',3));
%}

