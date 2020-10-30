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

WK= [1,.50,.75,.90,0];     %kk1
R=[.5,2,5,8,10];       %kk2

A1=[5;5;1]; A2=[5;3;1]; A3=[5;1;1];
AA= [A1,A2,A3];          %kk3
 
ANG= generateANG();
ANG= ANG*pi/180; %Converting angles from degrees to radians
%ANG= ANGG(:,1);%ang2,ang3];    %kk4
strain = [2,4,6,8];     %kk5
NNm=[1];                 %kk6 
NNn= [1];                %kk7
%*************************************************************************

%  generating 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 

 epsilonII=0.5;                                     % initial strain rate invariant  

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

             r = R(kk2);  % r here is the ratio of the inclusion viscosity to the matrix 
         
         
             
       
             %decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

              D = 0.5 * (L + L');
              W = 0.5 * (L - L');  
             
                         
               for kk3=1:(numel(AA)/3)
    
                % initial state of the ellipsoid
                                  
   
                        for kk4= 1:(numel(ANG)/3)
                          
                         
 
                            
                          
                            
                           
                           
                            for kk5=1:numel(strain)
                                
                                
                                for kk6=1:numel(NNm)
                                    
                                       
                                       for kk7=1:numel(NNn)     % this is the innermost loop !
                                           %  Power law coefficients for matrix and clast
                                           Nm=NNm(kk6);
                                           Nc=NNn(kk7);
                                           % calculations to follow:
                                           steps = 2000;   %  total steps of the computation
                                           strI=strain(kk5);
                                           %  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
                                           tincr = 0.05;
                                           % initial length of the ellipsoid axis
                                           a= AA(:,kk3);
                                           AAA=a; %noting initial shape;
                                           ang =   ANG(:,kk4);                      %Note: q and a are updated  in every calculation.(
                                           q = Q(ang);   % initial orientation matrix
                                           L_steps= zeros(3,3,steps);
                                           for k=1:steps
                                               d= q*D*q';              % Transforming D and W to Ellipsoidal axis system
                                               w= q*W*q'; 
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
                                               for k1=1:3
                                                   for l=1:3
                                                        S(k1,k1,l,l) = p(k1,k1)+ S_el(k1,k1,l,l);
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
                                               L_steps(:,:,k)= ll;
                                               
                                                                                              
                                               % update new a and q for next step
                                               a= qanew(:,4);
                                               q = qanew(1:3,1:3);
                                               [SI,Gamma,~]=Calc_Gamma(L,tincr,k);
                                               if SI>=strI
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
                                           %
                                           %----------------COMING BACK TO THE VPSC FOLDER----------------------------
                                           %---------------------------------------------------------------------------
                                            cd 'F:\Dropbox\VPSC_Modeling\Example\pureshear\';
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
                            
                                            fname = sprintf('LisoWk=%2.1f_r=%2.1f_Nm=%d_Nc=%d_shp_%d-%d-%d_IA=%d-%d-%d_gmma=%4.2f_SI=%4.2f.dat',Wk,r,Nm,Nc,AAA(1),AAA(2),AAA(3),ang(1),ang(2),ang(3),Gamma,SI);
                                            fileid= fopen(fname,'w');
                                            fprintf(fileid,formatspec1,steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' ');
                                            fprintf(fileid,formatspec2,'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', tincr);
                                            for iii=1:steps
                                                fprintf(fileid,formatspec,data(iii,:));
                                            end
                                            fclose(fileid);
                                            editvpsc7input(fname); %This function edits VPSC7.in file for VPSC input
                                                       
                                            command= 'F:\Dropbox\VPSC_Modeling\Example\pureshear\vpsc7.win64.exe';
                                            status= system(command);
                                            if Wk==0
                                                if r==.5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=0\r=.5\' fname];
                                                elseif r==2
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=0\r=2\' fname];
                                                elseif r==5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=0\r=5\' fname];
                                                elseif r==8
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=0\r=8\' fname];
                                                elseif r==10
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=0\r=10\' fname];
                                                end                                          
                                            elseif Wk==.5
                                                if r==.5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.5\r=.5\' fname];
                                                elseif r==2
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.5\r=2\' fname];
                                                elseif r==5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.5\r=5\' fname];
                                                elseif r==8
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.5\r=8\' fname];
                                                elseif r==10
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.5\r=10\' fname];
                                                end
                                            elseif Wk==.75
                                                if r==.5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.75\r=.5\' fname];
                                                elseif r==2
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.75\r=2\' fname];
                                                elseif r==5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.75\r=5\' fname];
                                                elseif r==8
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.75\r=8\' fname];
                                                elseif r==10
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.75\r=10\' fname];
                                                end
                                            elseif Wk==.90
                                                if r==.5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.90\r=.5\' fname];
                                                elseif r==2
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.90\r=2\' fname];
                                                elseif r==5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.90\r=5\' fname];
                                                elseif r==8
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.90\r=8\' fname];
                                                elseif r==10
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=.90\r=10\' fname];
                                                end
                                            elseif Wk==1
                                                if r==.5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=1\r=.5\' fname];
                                                elseif r==2
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=1\r=2\' fname];
                                                elseif r==5
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=1\r=5\' fname];
                                                elseif r==8
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=1\r=8\' fname];
                                                elseif r==10
                                                   fname1= ['F:\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=1\r=10\' fname];
                                                end                                 
                                            end
                                            % Plotting the CPO patterns using MTEX toolbox
                                            
                                            texfile= 'TEX_PH1.OUT';
                                            CS = crystalSymmetry('321',[4.9, 4.9, 5.4]);
                                            data = loadOrientation_generic(texfile,'CS',CS, 'ColumnNames', {'Euler1' 'Euler2' 'Euler3' 'Weights'},'Columns',[1,2,3,4],'Bunge');
                                            odf = calcODF(data);
                                            h=figure;
                                            plotPDF(odf,Miller(0,0,1,CS),'antipodal');
                                            setMTEXpref('xAxisDirection','East')
                                            set(gcf,'Visible','off');

                                            [path,name,~]= fileparts(fname1);
                                            fnamefig= [path '\' name '.fig'];
                                            savefig(fnamefig)
                                            
                                            %clear all variables in this
                                            %loop
                                            
                                       end
                                end
                                
                            end
                        end
               end
         end
   end
   
   movefile Liso* Iso_dat
