function Multi_Parameters_Aniso_HEM(Wk,r,a,ang,strI,m)

addpath(genpath('F:\Dropbox\MOPLA_Ankit_Matlab\'));
  
%  generating 4th-order identity tensors   
   [Jd, Js, Ja, ~, b] = Jnb(); 
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
                                

% for planar anisotropy parameter 
   Cm    = zeros(3,3,3,3);
     
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

         
      % r here is the ratio of the inclusion viscosity to the matrix Nn
 
 
         
   
        
        %decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

        D = 0.5 * (L + L');
        W = 0.5 * (L - L');  
      
   
        Ce=2*r*Jd;
        % Calculating Far field Stress Invariant for normalization :
        %Cm = C_bar_evl(:,:,:,:,200);      % stiffness tensor for the matrix
   
      
            
    
        % initial state of the ellipsoid
        % initial length of the ellipsoid axis
            
        AAA=a; %noting initial shape;
        Ang=ang*pi/180;
        q = Q(Ang);   % initial orientation matrix
        %  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
         tincr = 0.05;
        %  total steps of the computation
        steps = 1000;
       L_steps= zeros(3,3,steps);
        
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
         %
           if qa(1,4)/qa(3,4)>10
              qa(1,4)=0.5*qa(1,4);
           end

           if qa(2,4)/qa(3,4)>10
              qa(2,4)=0.5*qa(2,4);
           end
          qa = sortrows(qa,-4);  
                                           
          % updating L_ellipsoid
          ll= e+we;                          
          ll= q'*ll*q;                 
          L_steps(:,:,k)=ll;    
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
                            fname = sprintf('LanisoM=%dWk=%2.1fr=%2.1fshape%d-%d-%dthet=%d_phi=%d_gmma=%2.1f_SI=%4.2f.dat',m,Wk,r,AAA(1),AAA(2),AAA(3),ang(1),ang(2),Gamma,SI);
                            fileid= fopen(fname,'w');
                            fprintf(fileid,formatspec1,steps, 7, tincr, 298, 'nsteps', 'ictrl', 'eqincr' ,'temp',' ',' ',' ');
                            fprintf(fileid,formatspec2,'step', 'L11', 'L12', 'L13', 'L21', 'L22', 'L23', 'L31', 'L32', 'L33', tincr);
                            for iii=1:steps
                                 fprintf(fileid,formatspec,data(iii,:));
                            end
                            fclose(fileid);

                            
                            editvpsc7input(fname); %This function edits VPSC7.in file for VPSC input
                            
                            
                            command= 'F:\Dropbox\VPSC_Modeling\Example\pureshear\vpsc7.win64.exe';
                            system(command);
                                                       
                            if Wk==0
                                if r==.5
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=0\r=.5\' fname];
                                elseif r==2
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=0\r=2\' fname];
                                elseif r==5
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=0\r=5\' fname];
                                elseif r==8
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=0\r=8\' fname];
                                elseif r==10
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=0\r=10\' fname];
                                end                                          
                            elseif Wk==.5
                                if r==.5
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.50\r=.5\' fname];
                                elseif r==2
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.50\r=2\' fname];
                                elseif r==5
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.50\r=5\' fname];
                                elseif r==8
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.50\r=8\' fname];
                                elseif r==10
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.50\r=10\' fname];
                                end
                            elseif Wk==.75
                                if r==.5
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.75\r=.5\' fname];
                                elseif r==2
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.75\r=2\' fname];
                                elseif r==5
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.75\r=5\' fname];
                                elseif r==8
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.75\r=8\' fname];
                                elseif r==10
                                    fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.75\r=10\' fname];
                                end
                            elseif Wk==.90
                                if r==.5
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.90\r=.5\' fname];
                                elseif r==2
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.90\r=2\' fname];
                                elseif r==5
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.90\r=5\' fname];
                                elseif r==8
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.90\r=8\' fname];
                                elseif r==10
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=.90\r=10\' fname];
                                end
                             elseif Wk==1
                                if r==.5
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=1\r=.5\' fname];
                                elseif r==2
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=1\r=2\' fname];
                                elseif r==5
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=1\r=5\' fname];
                                elseif r==8
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=1\r=8\' fname];
                                elseif r==10
                                   fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\AnisotropicHEM\Wk=1\r=10\' fname];
                                end                                 
                            end
                                
                             % Plotting the CPO patterns using MTEX toolbox
                                            
                                            texfile= 'TEX_PH1.OUT';
                                            CS = crystalSymmetry('321',[4.9, 4.9, 5.4]);
                                            data = loadOrientation_generic(texfile,'CS',CS, 'ColumnNames', {'Euler1' 'Euler2' 'Euler3' 'Weights'},'Columns',[1,2,3,4],'Bunge');
                                            odf = calcODF(data);
                                            figure;
                                            plotPDF(odf,Miller(0,0,1,CS),'antipodal');
                                            setMTEXpref('xAxisDirection','East')
                                            %set(gcf,'Visible','off');
                                            figure;
                                            plotPDF(odf,Miller({1,0,-1,0},{0,1,-1,1},{2,-1,-1,0},{1,0,-1,1},CS),'antipodal')
                                            setMTEXpref('xAxisDirection','East')
                                           

                                            %[path,name,~]= fileparts(fname1);
                                            %fnamefig= [path '\' name '.fig'];
                                           % savefig(gcf,fnamefig,'compact')
                                           % filenamenew= [path '\' name, '.jpg'];
                                           % saveas(gcf,filenamenew)
                                           % close(gcf)
                                            movefile Laniso* Aniso_dat                                        
                                       
end
   
  








