% Generating variable velocity gradient tensors for different initial parameters
% Simulating motion of single ellipsoid in a planar  anisotropic HEM 
 

% Different initial parameters

% Wk (Kinematic vorticity number) ranging from 0-1. Wk= 0,.5,.75,.90,1
% Viscosity ratio of RDE to HEM's Nn (since HEM is planar anisotropic) r=.5, 1.5, 2,5,10
% Inital orientation of the RDE given by theta 0,30,60,90
% Initial shape of the RDE  10:10:1, 10:5:1, 10:2:2(Oblate, triaxial and Prolate)

%*************************************************************************
% NOTE: these are the following parameters that are changed in the nested
% for-loops to run simulations for various initial conditions

MM= [1];%,2,5,10]; %kk0
WK= [0];     %kk1
R=[5];%1.5,2,5,10];       %kk2
A1=[10;10;1]; A2=[10;5;1]; A3=[10;2;1.99];
AA= [A2];%,A2,A3];          %kk3
ttheta= [30];%,30,60];    %kk4
strain = [4];%,4,6,8];
%*************************************************************************

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
   
for kk0=1: numel(MM)
   % for planar anisotropy parameter 
   Cm    = zeros(3,3,3,3);
         m= MM(kk0);
         Nn=1;

         for i=1:3
             for j=1:3
                  for kt=1:3
                      for l=1:3
                           if rem(i+j,2)==0
                              Cm(i,j,kt,l)= 2*Nn*Jd(i,j,kt,l);
                           else
                              Cm(i,j,kt,l)= 2*Nn*Jd(i,j,kt,l)/m;
                           end
                      end
                  end
             end
         end 
   for kk1=1:numel(WK)
    

   

         %  Input parameters--------------------------------------------------------
         Wk= WK(kk1);
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

         for kk2=1:numel(R)
         % defining the stiffness tensor of the planar anisotropic HEM 

         r = R(kk2);  % r here is the ratio of the inclusion viscosity to the matrix Nn
 
 
         
   
         %load('AnisotropicC_200inclusionsB','C_bar_evl');
        %decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

        D = 0.5 * (L + L');
        W = 0.5 * (L - L');  
        Wvec= [2*W(3,2); 2*W(1,3);2*W(2,1)];  % vorticity vector
        zaxis= [0;0;1];
        dotprod1= dot(Wvec,zaxis);

   
        Ce=2*r*Jd;
        % Calculating Far field Stress Invariant for normalization :
        %Cm = C_bar_evl(:,:,:,:,200);      % stiffness tensor for the matrix
   
        Sigma = Multiply(Cm,D);        % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
        Sigma_Inva = inva(Sigma);
        strainInva = inva(D);
               for kk3=1:(numel(AA)/3)
    
                % initial state of the ellipsoid
                                     % initial length of the ellipsoid axis
                a= AA(:,kk3);
                AAA=a; %noting initial shape;
   
                        for kk4= 1:numel(ttheta)
                            ii= ttheta(kk4); % i is the initial orientation angle
                            ang = [0; 0;(pi/180)*ii];                        % See Jiang 2007a for choice of angles.
    
                            q = Q(ang);   % initial orientation matrix
 
                            %  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
                            tincr = 0.05;
                            %  total steps of the computation
                            steps = 1000;
          
                            % calculations to follow:
                            dotprodcheck= zeros(steps,1);  %checking the sense of vorticity  
 
                            L_steps= zeros(3,3,steps);
                            
                            for kk5=1:numel(strain)
                                steps = 1000;
                                strI=strain(kk5);
                                for k=1:steps
                                    
   
                                    d= q*D*q';                                    % Transforming D and W to Ellipsoidal axis system
                                    w= q*W*q'; 
                                    Cm_a    = Transform(Cm, q);
                                    Carray  = C2OneDarray(Cm_a);
                                    T= TGreen(a,  Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                                                   Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
                                    S = Contract(Jd, Contract(T,Cm_a));
                                    PI= Contract(Ja, Contract(T,Cm_a));
   
                                    invS = Inverse((S),b);
                                    h= Contract(Cm_a, invS);                        
    
                                    H = Contract(Cm_a,(invS-Js));           % Eq.12a Jiang 2014
   
                                    A = Contract(FourTensorInv(H + Ce),(H + Cm_a));                     
                                    e= Multiply(A,d);    % e is strain rate in ellipsoid
                                    u1= e-d;                  
                                                                    
                                   %  vorticity of RDE        
       
                                    uu1    = Contract(PI,invS);
                                    we    = Multiply(uu1, u1)+ w; 
                                    wEp   = Wd(a, w, e);
             
                                    %  angular velocity of RDE
                                    theta = we - wEp; 
                                    %  update Q
                                    qq   = (RodrgRot(-theta * tincr)) * q; 
                                    %  update a
                                    aa   = a.* exp(diag(e) * tincr); 
                                    %  make sure that Q and a are in the descending oreder of a(a1>=a2>=a3)          
                                    qa    = sortrows([qq aa],-4);          
                                    %  Boudinage if the RDE is too elongated or flattened(a1:a3>100 or a2:a3>100)
                                      %{
                                       if qa(1,4)/qa(3,4)>10
                                          qa(1,4)=0.5*qa(1,4);
                                       end

                                       if qa(2,4)/qa(3,4)>10
                                          qa(2,4)=0.5*qa(2,4);
                                       end
                                           qa = sortrows(qa,-4);  
                                           %}
                                    % updating L_ellipsoid
    
                                    ll= e+we;
                                    ll= q'*ll*q;
                                    L_steps(:,:,k)=ll;
                                    wEE= q'*we*q;
                                    wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
                                    dotprod= dot(wvec,zaxis);
                                    dotprodcheck(k,1)= dotprod; % checking the shear sense !
        
        
                                    % updated for the next step     
                                    a   = qa(:,4);
                                    q = qa(1:3,1:3);   
                                    [SI,Gamma,~]=Calc_Gamma(L,tincr,k);
                                    if SI>=strI   % This line decides the strain intensity of the macroscale deformation
                                       break
                                    end
                                end
                            steps=k;

                            numb= ones(steps,1);
                            tcr= tincr*ones(steps,1);
                            Lstep= zeros(steps,9);

                            for k1=1:steps
                                numb(k1)=k1;
                                kk= L_steps(:,:,k1);
                                Lstep(k1,:)= reshape(kk',1,9);
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
                            
                            % temporary work around the bug
                            if Wk==0
                                Gamma=1;
                            end
                            
                            fname = sprintf('L_aniso_m=%d_Wk=%2.1f_r=%d_shape_%d-%d-%d_theta=%d_gamma=%4.2f_SI=%4.2f.dat',MM(kk0),Wk,r,AAA(1),AAA(2),AAA(3),ii,Gamma,SI);
                            fileid= fopen(fname,'w');
                            fprintf(fileid,formatspec1,steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' ');
                            fprintf(fileid,formatspec2,'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', tincr);
                            for iii=1:steps
                                 fprintf(fileid,formatspec,data(iii,:));
                            end
                            fclose(fileid);

                            
                            editvpsc7input(fname); %This function edits VPSC7.in file for VPSC input
                            
                            
                            command= 'D:\Ankit\Dropbox\VPSC_Modeling\Example\pureshear\vpsc7.win64.exe';
                            status= system(command);

                            % Plotting the CPO patterns using MTEX toolbox
                            PlotTexture_AB(fname);
                            %{
                            fname = sprintf('test_L_aniso_Wk=%d_r=%d_shape_%d-%d-%d_theta=%d.xls',Wk,r,AA(1,kk3),AA(2,kk3),AA(3,kk3),ii);
                            firstlin= {steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' '};
                            seclin= {'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', 'tincr'};
                            datanew= {firstlin;seclin;data};
                            xlswrite(fname,firstlin,'sheet1');
                            xlswrite(fname,seclin,'sheet1',sprintf('A%d',2));
                            xlswrite(fname,data,'sheet1',sprintf('A%d',3));
                            %}
                            end
                        end
               end
         end
   end
   
end
