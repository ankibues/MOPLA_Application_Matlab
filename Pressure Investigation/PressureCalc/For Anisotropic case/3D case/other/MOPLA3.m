%function MOPLA()
% same code for running in different system
%  MOPLA.m        update:May 20,2016
%
% Modeling the mechanicall behavior of the heterogeneous rock mass which is
% composed of rheologically distinctive elements(RDEs). Here, the heterogeneous
% rock mass is replaced by a homogeneous effective medium(HEM) and the
% rehological properties of HEM at a point are represented by the overall
% effective propertives of all surrounding RDEs.
%
% Scripts are based on the self-consistent solution of the partitioning and
% homogenization equations.(Jiang, 2013; Jiang, 2014)
% 
%  clear all variables, Comment Window and figures-------------------------  
   
   
%  Input parameters--------------------------------------------------------
%  Imposed macroscale flow filed
   L     = [0 1 0 ;0 0 0;0 0 0]; 
%  the number of RDEs(Rheologically Distinctive Elements) 
   n     = 10;
%  the maximum semi-axils of RDEs (a1>a2>a3; a3=1)
   a1    = 10;
%  the maximum viscosity of RDEs
   vmax  = 10;
%  the minimum viscosity of RDEs
   vmin  = 1;
%  the maximum stress exponent of RDEs
   Nmax  = 4;
%  the minimum stress exponent of RDEs
   Nmin  = 1;
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;
%  total steps of the computation
   steps = 20;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');   
%  generate 4th-order identity tensors   
   [Jd, ~, Ja, ~] = FourIdentity();
   
%  obtain weights and nodes before the loop // these won't be needed for
%  isotropic case
   gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
   
%  generate a population of uniformly distributed RDEs with
% a1:a2:1,a2 from a1 to 1.

  [a, ang] = RandAANG(a1,n);
%  generate the viscosities(eta) of RDEs
  % eta = vmin + (vmax-vmin)*rand(1,n);
    eta = 5*ones(1,n);                                
%  generate the stress exponents(Ne) of RDEs
    Ne  = Nmin + (Nmax-Nmin)*rand(1,n);
   
%  generate the transfermation tensors(q),compliance                                                                                           tensors(m)and pre-strain
%  rate term(e0) of RDEs
   q   = zeros(3,3,n);
   m   = zeros(3,3,3,3,n);
   e0  = zeros(3,3,n);
   for i= 1:n
        
       q(:,:,i)= Q(ang(:,i));
       m(:,:,:,:,i) = 1/(2*eta(1,i))*Jd;
       e0(:,:,i) = (1-Ne(1,i))*D;
   end
%  strain rate invariant at which ellipsoid viscosity is defiend  
   REF  = Inva(D)*ones(1,n);                           
   
%  initial guess of the properties of the HEM:-----------------------------
   bulkstresscheck= zeros(1,steps);
   bulkstrainratecheck=zeros(1,steps);
%  generate the tangent compliances(mtan) of RDEs
   mtan   = zeros(3,3,3,3,n);
   for i=1:n
       mtan(:,:,:,:,i)= Ne(1,i)* m(:,:,:,:,i);         
   end
%  initial homogenized compliance(M_bar) of HEM
   M_bar            = 1/n*sum(mtan,5);                 
%  initial homogenized stiffness(C_bar) of HEM
   C_bar            = FourTensorInv(M_bar);
%  initial pre-strain rate term(E0) of HEM   
   mm               = FourTensorInv(1/n*sum(m,5));     
   t1               = Contract(M_bar,mm);              
   E0               = Multiply((Jd-t1),D); 
%  the macroscopic stress(Stress_bar) of HEM
   Stress_bar       = Multiply(C_bar,D-E0);     
%  far-field stress 
   SIGMA            = Stress_bar;        
%  invariant of the HEM deviatoric stress(macroscale)             
   BStress          = Inva(Multiply(Jd,Stress_bar));   
%  partitioned stress fields inside RDEs
   sigma            = repmat(Stress_bar,1,1,n);       
   
%  preallocate variables:--------------------------------------------------
        c_bar = zeros(3,3,3,3,n); 
        T     = zeros(3,3,3,3,n); 
        H_arc = zeros(3,3,3,3,n);
        e     = zeros(3,3,n);     
        eI    = zeros(1,n);        
        vis   = zeros(1,n);         
        Carray = zeros(n,81);
        Q_evl     = zeros(3,3,n,steps);
        A_evl     = zeros(3,n,steps);
        C_bar_evl = zeros(3,3,3,3,steps);
        M_bar_evl = zeros(3,3,3,3,steps);
        stress_evl = zeros(3,3,steps);
        stress_evl2 = zeros(3,3,steps);
 %  A Self-consistent approach---------------------------------------------
 for l=1:steps   
    tic
     for i=1:100
     
            %   Outer loop:When all elements are calculated, using the new set
            %   of M(k),B(k),e0(k),beta(k) update new M_bar and E0 until it
            %   approaches balance.
        
        
        MB = zeros(3,3,3,3);  % for <M:B>
        BB = zeros(3,3,3,3);  % for <B>
        Mb = zeros(3,3);      % for <M:B+e0>
        bb = zeros(3,3);      % for <beta>
%  describe C in the RDE's coordinate system 
       c_bar = Transform2(C_bar,q,n);       % check how much time is improved using this, in an actual run !
       
        for r=1:n    %this loop  should be parallelized !
            %   S and H_arc for each RDP are calculated starting with
            %   the M_bar and E0 defined at the homogeneous macrosclae
            %   strain-rate state.
            

%  rewrite the matrix stiffness tensor into a 1D array format 
            Carray(r,:)  = C2OneDarray(c_bar(:,:,:,:,r));   
%  compute the 4th-order Green tensor(T) (Vectorize this/ Done ! )            
            T(:,:,:,:,r)     = TGreen(a(:,r),  Carray(r,:), Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
%  the 4th-order Interaction tensor(H_arc) (vectorize this)
            s1               = Contract(Jd, Contract(T(:,:,:,:,r),c_bar(:,:,:,:,r)));
            s                = Transform(s1,q(:,:,r)'); %% is this really needed ?? Yes, this transforms back to the general coordinate system.                          
            H_arc(:,:,:,:,r) = Contract(FourTensorInv(FourTensorInv(s)-Jd),M_bar);  
%  record the initial partitioned stress of RDE at the current homogeneous macroscopic state             
            stress           = sigma(:,:,r);
           
            for k=1:100
                %   Inner loop: B and beta(b) for each RDP are calculated starting with
                %   the M_bar and E0 defined at the homogeneous macrosclae
                %   strain-rate state uitil is appraches the balance.

%  the 4th-order stress-partitioning tensor(B)                
                t0              = Ne(1,r)*m(:,:,:,:,r);                            
                t1              = FourTensorInv((t0+H_arc(:,:,:,:,r)));                   
                B               = Contract(t1,(M_bar+H_arc(:,:,:,:,r)));           
                t2              = E0-e0(:,:,r);
%  the second order stress-partitioning tensor(beta)
                beta            = Multiply(t1,t2);   
%  calculate new partitioned stress field inside RDE
                sigma(:,:,r)    = Multiply(B,SIGMA)+beta;                          
%  compare the current partitioned stress and the previous one                 
                alpha           = abs(Inva(sigma(:,:,r)-stress)/Inva(stress));
%  replace the previous partitioned stress with the current one
                stress          = sigma(:,:,r);
%  calculate new partitioned strain rate of RDE               
                e(:,:,r)        = Multiply(m(:,:,:,:,r),sigma(:,:,r));                    
%  calculate new pre-strain rate term of RDE                
                e0(:,:,r)       = (1-Ne(1,r))*e(:,:,r);
%  calculate new strain rate invariant
                e1              = Inva(e(:,:,r));
%  calculate new viscosity of RDE at the new strain rate state               
                t3              = 1/Ne(1,r)-1;
                v1              = (e1/REF(1,r))^t3 *eta(1,r);     
%  calculate new compliance tensor of RDE
                m(:,:,:,:,r)    = 1/(2*v1)*Jd;                                     
               
%  The inner loop for an RDE terminates when the current stress coincides
%  with the previous one within a specific tolerance
               if alpha<0.01 
                   break
               end   
               
            end
%  update the strain-rate invariant for an RDE           
            eI(1,r)  = e1;  
%  update the effective viscosity for an RDE
            vis(1,r) = v1; 
            
            MB = MB+Contract(t0,B);                    %%<M(tan):B>*n%%
            Mb = Mb+Multiply(t0,beta)+ e0(:,:,r);       %%<M(tan):beta+e0>*n%%
            BB = BB+B;                                 %%<B>*n$$
            bb = bb+beta;                              %%<beta>*n%%
        end
%  calculate new homogenized compliance for HEM         
        invB         = FourTensorInv(BB/n);
        New_M_bar    = Contract(MB/n,invB);  
%  compare the current homogenized compliance and the previous one        
        delta        = Norm(New_M_bar-M_bar)/Norm(M_bar);
%  replace the previous homogenized compliance with the current one        
        M_bar        = New_M_bar;  
%  calculate new homogenized stiffness for HEM
        C_bar        = FourTensorInv(M_bar);  
%  calculate new pre-strain rate term for HEM
        E0           = (Mb/n)-Multiply(M_bar,(bb/n));     
%  calculate new homogeneous macroscopic stress for HEM
        Stress_bar   = Multiply(C_bar,D-E0);  % ********************   E3a (jiang 2014)         
%  calculate new far-filed stress
        SIGMA        = Multiply(invB,(Stress_bar- bb/n)); % looks like average of E5 to me? discuss !
%  calculate new secound invariant of the macroscopic deviatoric stress  
        NewStress    = Inva(Multiply(Jd,Stress_bar));  
%  compare the current macroscopic deviatoric stress and the previous one
        gamma        = abs(NewStress/BStress-1);
%  replace the previous invariant of the macroscopic deviatoric stress with the current one         
        BStress      = NewStress; 
        
%  The outer loop continues until the current homogenized compliance and 
%  macroscopic deviatoric stress coincide with the previous ones respectively
%  within specific tolerances        
        if delta<0.02 && gamma<0.02
            break
       end
       
     end
     bulkstresscheck(:,l)= NewStress; % checking how bulk strain rate and stress evolves !
     stress_evl(:,:,l)= Stress_bar;
     stress_evl2(:,:,l)= SIGMA;
%  update the strain-rate invar.iants for all RDEs
     REF=eI;    
%  update the effective viscosities for all RDEs
     eta=vis;   
%  write updated C_bar to C_bar_evl
     C_bar_evl(:,:,:,:,l)=C_bar;
%  write updated M_bar to M_bar_evl
     M_bar_evl(:,:,:,:,l)=M_bar;
%  Evolution of RDEs-------------------------------------------------------   
     for j=1:n
%  describe D,W in the RDE's coordinate system 
        D_bar = q(:,:,j)*D*q(:,:,j)'; 
        W_bar = q(:,:,j)*W*q(:,:,j)'; 
%  calculate Eshelby Tensors(S, PI) based on T
        z    = Contract(T(:,:,:,:,j),c_bar(:,:,:,:,j));        
        S    = Contract(Jd,z);                                
        PI   = Contract(Ja,z);        
%  strain rate of RDE
        de    = q(:,:,j)*e(:,:,j)*q(:,:,j)'; %Ellipsoid strain-rate
%  vorticity of RDE        
        invS  = FourTensorInv(S);
        u1    = Contract(PI,invS);
        u2    = de - D_bar;
        we    = Multiply(u1, u2)+ W_bar; 
        wEp   = Wd(a(:,j), W_bar, de);        
%  angular velocity of RDE
        theta = we - wEp; 
%  update Q
         qq   = (RodrgRot(-theta * tincr)) * q(:,:,j); 
%  update a
         aa   = a(:,j).* exp(diag(de) * tincr); 
%  make sure that Q and a are in the descending oreder of a(a1>=a2>=a3)          
         qa    = sortrows([qq aa],-4);          
%  Boudinage if the RDE is too elongated or flattened(a1:a3>100 or
%  a2:a3>100) changing this critical length
         if qa(1,4)/qa(3,4)>25
            qa(1,4)=0.5*qa(1,4);
        end

        if qa(2,4)/qa(3,4)>25
            qa(2,4)=0.5*qa(2,4);
        end
        qa = sortrows(qa,-4);       
        a(:,j)   = qa(:,4);
        q(:,:,j) = qa(1:3,1:3);
%  write updated Q to Q_evl        
        Q_evl(:,:,j,l)=q(:,:,j);
%  write updated a to A_evl
        A_evl(:,j,l)=a(:,j);
        
     end
  toc
  disp('no of Steps completed');
  disp(l);  
  
 end 
   
  % save Q_evl, A_evl, C_bar_evl and M_bar_evl to the current workspace   
 %save('AnisotropicC_100inclusions.mat','Q_evl','A_evl','C_bar_evl','M_bar_evl','stress_evl','stress_evl2','bulkstresscheck');
 
%  get the final orientaions of the rigid ellipsoids
   Q_final = Q_evl(:,:,:,steps);

%  Equal-area projection

%     compute two spherical angles for three axes
      [a1_ang, a2_ang, a3_ang] = ConvertQ2Angs(Q_final);

      
%     compute r for equal-area projection, both hemispheres will be plotted
%     a1
      [~,a1in]  = find(a1_ang(2,:)<=(0.5*pi));
      [~,a1out] = find(a1_ang(2,:)>(0.5*pi));
      r1(a1in)  = sqrt(2) * sin(a1_ang(2,a1in)./2);
      r1(a1out) = sqrt(2) * cos(a1_ang(2,a1out)./2);
%     a2
      [~,a2in]  = find(a2_ang(2,:)<=(0.5*pi));
      [~,a2out] = find(a2_ang(2,:)>(0.5*pi));
      r2(a2in)  = sqrt(2) * sin(a2_ang(2,a2in)./2);
      r2(a2out) = sqrt(2) * cos(a2_ang(2,a2out)./2);
%     a3       
      [~,a3in]  = find(a3_ang(2,:)<=(0.5*pi));
      [~,a3out] = find(a3_ang(2,:)>(0.5*pi));
      r3(a3in)  = sqrt(2) * sin(a3_ang(2,a3in)./2);
      r3(a3out) = sqrt(2) * cos(a3_ang(2,a3out)./2);
   
%     equal-area projections of a1, a2, a3   
%     a1      
      subplot(1,3,1);
      t = 0 : .01 : 2 * pi;
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
%     phi<=pi/2, plot red dots 
      polar(a1_ang(1,a1in),r1(a1in),'.r')
%     phi>pi/2, plot red dots 
      polar(a1_ang(1,a1out),r1(a1out),'.r')
      hold off
      title('a1')
      
%     a2   
      subplot(1,3,2);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a2_ang(1,a2in),r2(a2in),'.r')
      polar(a2_ang(1,a2out),r2(a2out),'.r')
      hold off
      title('a2')
      
%     a3   
      subplot(1,3,3);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a3_ang(1,a3in),r3(a3in),'.r')
      polar(a3_ang(1,a3out),r3(a3out),'.r')
      hold off
      title('a3')
   
% Flinn diagram
    A = A_evl(:,:,steps);
    x = log(A(2,:)./A(3,:));
    y = log(A(1,:)./A(2,:));
    
    figure('Name','Flinn diagram: The shape of ellipsoids');
    
    plot(x,y,'.r',0:5,0:5,'-k')
    xlabel('ln(a2/a3)')
    ylabel('ln(a1/a2)')
    axis square
 % plotting bulk stress and strain rate evolution
 bulkstressstep1= bulkstresscheck(:,1);
 bulkstresscheck= bulkstresscheck/bulkstressstep1;
 X=1:steps;
 figure
 plot(X,bulkstresscheck)
 
  save('AnisotropicC_100inclusions.mat','Q_evl')
    
%end
  
   
 