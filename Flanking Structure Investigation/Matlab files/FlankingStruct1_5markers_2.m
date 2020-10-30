% Deformation of marker with exterior points around inclusion- modelling a
% flanking structure

% Note : In this code, mesh is defined with respect to external coordinate system

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
 
%  Input parameters--------------------------------------------------------
%  Matrix flow field
  % L     = [0 1 0;0 0 0;0 0 0]; 
D = 0.5 * (L + L');
% normalizing with respect to strain rate
L= L/norm(D);

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps =50;
   
%  Power law coefficients for matrix and inclusion
   Nm=1;
   Nc=1; 
   r= .0001;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% Stiffness tensor of the matrix   
   Cm    = 2*Jd;
   
% initial state of the ellipsoid                                    
   a= [1000;.5;.1];
   A=a;% initial length of the ellipsoid axis
   pp= 3*pi/4;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
 
% Here, I define the markers in the external coordinate system

% 3D meshgrid in external coordinate
        xgv     = -3:.1:3;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     =0*ones(size(xgv));           % grid vector: y'axis,a2
        zgv     = zeros(size(xgv)); %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        ep    = cat(1,xgv,ygv,zgv);
        [~,N]= size(ep);
       %{
        xgv_1     = -3:.03:3;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv_1     = 0*ones(size(xgv));           % grid vector: y'axis,a2
        zgv_1     = zeros(size(xgv)); %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        ep_1    = cat(1,xgv_1,ygv_1,zgv_1);
        
        
        xgv_2     = -3:.03:3;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv_2     = -.2*ones(size(xgv));           % grid vector: y'axis,a2
        zgv_2     = zeros(size(xgv)); %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        ep_2    = cat(1,xgv_2,ygv_2,zgv_2);
 
        
        xgv_3     = -3:.03:3;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv_3     = .6*ones(size(xgv));           % grid vector: y'axis,a2
        zgv_3     = zeros(size(xgv)); %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        ep_3    = cat(1,xgv_3,ygv_3,zgv_3);
        
        
        xgv_4     = -3:.03:3;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv_4     = -.6*ones(size(xgv));           % grid vector: y'axis,a2
        zgv_4     = zeros(size(xgv)); %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        ep_4    = cat(1,xgv_4,ygv_4,zgv_4);
        
        %}
        
        
        
        
      % bringing markers into interior coordinate system
    
    parfor rrr=1:N                               
        ep(:,rrr)= q*ep(:,rrr);
        %{
        ep_1(:,rrr)= q*ep_1(:,rrr);
        ep_2(:,rrr)= q*ep_2(:,rrr);
        ep_3(:,rrr)= q*ep_3(:,rrr);
        ep_4(:,rrr)= q*ep_4(:,rrr);
        %}
    end
    
    

    
% frames to capture in every step:
  
  F(steps)= struct('cdata',[],'colormap',[]);
  
 %

for i=1:steps
   tic
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
    
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
    
    
    invS= Inverse(S,b);
    h= Contract(2*Jd/Nm, invS);           % Here is the Cinverse, divided by Nm             
    p_in = R_Multiply(p,h);
    E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
    r_new= E(1,4);
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    ll= dE + we;                                    % Partitioned velocity gradient tensor


%
    % Also, I want to check if the marker points are outside the ellipsoid.
    
    %------------------------for first marker
    indd= (ep(1,:)/a(1)).^2 + (ep(2,:)/a(2)).^2 + (ep(3,:)/a(3)).^2>1;
    ep= ep(:,indd);     % this will shortlist all the exterior points of the marker.
     [~,N]= size(ep);


     %{
    %------------------------for second marker
    indd1= (ep_1(1,:)/a(1)).^2 + (ep_1(2,:)/a(2)).^2 + (ep_1(3,:)/a(3)).^2>1;
    ep_1= ep_1(:,indd1);
     [~,N1]= size(ep_1);
    %------------------------for third marker
    indd2= (ep_2(1,:)/a(1)).^2 + (ep_2(2,:)/a(2)).^2 + (ep_2(3,:)/a(3)).^2>1;
    ep_2= ep_2(:,indd2);  
     [~,N2]= size(ep_2);
    %------------------------for fourth marker
    indd3= (ep_3(1,:)/a(1)).^2 + (ep_3(2,:)/a(2)).^2 + (ep_3(3,:)/a(3)).^2>1;
    ep_3= ep_3(:,indd3);  
     [~,N3]= size(ep_3);
     %------------------------for fifth marker
    indd4= (ep_4(1,:)/a(1)).^2 + (ep_4(2,:)/a(2)).^2 + (ep_4(3,:)/a(3)).^2>1;
    ep_4= ep_4(:,indd4);  
     [~,N4]= size(ep_4);
     
     
     %}
     
     
%}
  
 
    % Now, once all the exterior points are available , we calculate
    % velocity field and displacement field for them.
   %------------------------------------------------------ 
    % for first marker !
    vel_f= zeros(2,1,N);
    xxx= ep(2,:);
    yyy= ep(3,:);
    xxx(isnan(xxx)==1)=[];
    yyy(isnan(yyy)==1)=[];
    % [~,N]=size(ep);
    %
    parfor ir=1:N
        vel_f(:,:,ir)= Velocity_field_calc4(xxx(ir),yyy(ir),a,invS,Jd,d,u2,w,ll);
       
    end
    if i>(steps/2)
        tincr=-tincr;
    end
    disp_f= vel_f*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form
    
    kkk= zeros(1,N);
    DISP_F= [kkk; squeeze(disp_f)];
    newpos_plot= zeros(size(DISP_F));
    newpos= ep + DISP_F;  % this is the new position calculated for the marker points.

    %{
    %{ note : this is still calculated in inclusion coordinate system. Has to be converted to external coordinate system. This
    % is done later, after updating the q for the next step.
    %
  %------------------------------------------------------ 
    % for second marker !
   vel_f1= zeros(2,1,N1);
    xxx1= ep_1(2,:);
    yyy1= ep_1(3,:);
    xxx1(isnan(xxx1)==1)=[];
    yyy1(isnan(yyy1)==1)=[];
     [~,N1]=size(ep_1);
    parfor ir=1:N1
        vel_f1(:,:,ir)= Velocity_field_calc4(xxx1(ir),yyy1(ir),a,invS,Jd,d,u2,w,ll);
       
    end
    disp_f1= vel_f1*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form  
    kkk1= zeros(1,N1);
    DISP_F1= [kkk1; squeeze(disp_f1)];
    newpos_plot1= zeros(size(DISP_F1));
    newpos1= ep_1 + DISP_F1;  % this is the new position calculated for the marker points.
                   
         
     
   %------------------------------------------------------ 
    % for third marker !
    vel_f2= zeros(2,1,N2);
    xxx2= ep_2(2,:);
    yyy2= ep_2(3,:);
    xxx2(isnan(xxx2)==1)=[];
    yyy2(isnan(yyy2)==1)=[];
     [~,N2]=size(ep_2);
    parfor ir=1:N2
        vel_f2(:,:,ir)= Velocity_field_calc4(xxx2(ir),yyy2(ir),a,invS,Jd,d,u2,w,ll);
       
    end
    disp_f2= vel_f2*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form  
    kkk2= zeros(1,N2);
    DISP_F2= [kkk2; squeeze(disp_f2)];
    newpos_plot2= zeros(size(DISP_F2));
    newpos2= ep_2 + DISP_F2;  % this is the new position calculated for the marker points.
   %--------------------------------------   
   
   %------------------------------------------------------ 
    % for fourth marker !
    vel_f3= zeros(2,1,N3);
    xxx3= ep_3(2,:);
    yyy3= ep_3(3,:);
    xxx3(isnan(xxx3)==1)=[];
    yyy3(isnan(yyy3)==1)=[];
     [~,N3]=size(ep_3);
    parfor ir=1:N3
        vel_f3(:,:,ir)= Velocity_field_calc4(xxx3(ir),yyy3(ir),a,invS,Jd,d,u2,w,ll); 
    end
    disp_f3= vel_f3*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form  
    kkk3= zeros(1,N3);
    DISP_F3= [kkk3; squeeze(disp_f3)];
    newpos_plot3= zeros(size(DISP_F3));
    newpos3= ep_3 + DISP_F3;  % this is the new position calculated for the marker points.
   %--------------------------------------
   
   %------------------------------------------------------ 
    % for fifth marker !
    vel_f4= zeros(2,1,N4);
    xxx4= ep_4(2,:);
    yyy4= ep_4(3,:);
    xxx4(isnan(xxx4)==1)=[];
    yyy4(isnan(yyy4)==1)=[];
     [~,N4]=size(ep_4);
    parfor ir=1:N4
        vel_f4(:,:,ir)= Velocity_field_calc4(xxx4(ir),yyy4(ir),a,invS,Jd,d,u2,w,ll);
       
    end
    disp_f4= vel_f4*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form  
    kkk4= zeros(1,N4);
    DISP_F4= [kkk4; squeeze(disp_f4)];
    newpos_plot4= zeros(size(DISP_F4));
    newpos4= ep_4 + DISP_F4;  % this is the new position calculated for the marker points.
   %--------------------------------------
   
   %}
     
       
     parfor rrr=1:N                                % bringing to general coordinate system
         newpos_plot(:,rrr)= q'*newpos(:,rrr);
     end
         
    %----------------------------
     %{
     parfor rrr=1:N1                                % bringing to general coordinate system
         newpos_plot1(:,rrr)= q'*newpos1(:,rrr);
     end   
     
     parfor rrr=1:N2                                % bringing to general coordinate system
         newpos_plot2(:,rrr)= q'*newpos2(:,rrr);
     end  
     
      parfor rrr=1:N3                                % bringing to general coordinate system
         newpos_plot3(:,rrr)= q'*newpos3(:,rrr);
      end  
      
      parfor rrr=1:N4                                % bringing to general coordinate system
         newpos_plot4(:,rrr)= q'*newpos4(:,rrr);
      end  
     
     %}
     
     
     
     %------------------------------------------------------------
 
    
           
        % 3D meshgrid in internal coordinate -----interior points of
        % inclusion-------------------------------------------------------------------------
        
        xgv1     = 0; %-5:.05:5;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv1    = -2:.1:2;           % grid vector: y'axis,a2
        zgv1     = -2:.1:2;  % grid vector: z'axis,a3
        [X1,Y1,Z1] = meshgrid(xgv1,ygv1,zgv1);
        ind1= (X1./a(1)).^2 + (Y1./a(2)).^2 + (Z1./a(3)).^2;
        x_ex1  = X1(ind1<1);
        y_ex1  = Y1(ind1<1);
        z_ex1  = Z1(ind1<1);
        epp    = cat(1,x_ex1',y_ex1',z_ex1'); % These are points corresponding the inclusion
        [~,NN]= size(epp);
        epp1= zeros(size(epp));
        for rrr=1:NN                                % bringing to general coordinate system
           epp1(:,rrr)= q'*epp(:,rrr);    % in exterior coordinate system for further conversion
        end

        qanew= QL(a,q,theta,dE,tincr);  
         a=qanew(:,4);
         q=qanew(1:3,1:3);
        
       % plotting part for each step
   
       XX= newpos_plot(1,:);
       YY= newpos_plot(2,:); % first marker 
       %{
       XX1= newpos_plot1(1,:);
       YY1= newpos_plot1(2,:); % second marker 
       XX2= newpos_plot2(1,:);
       YY2= newpos_plot2(2,:); % third marker 
        XX3= newpos_plot3(1,:);
       YY3= newpos_plot3(2,:); % fourth marker 
        XX4= newpos_plot4(1,:);
       YY4= newpos_plot4(2,:); % fifth marker 
       %}
       
      
       XXX= epp1(1,:);
       YYY= epp1(2,:);  % inclusion
        f= figure('visible','off');
      
        plot(XX,YY,'.','MarkerSize',2)     % this is for the marker 1 !
        hold on
       % plot(XX1,YY1,'.','MarkerSize',2)     % this is for the marker 2 !
        hold on
       % plot(XX2,YY2,'.','MarkerSize',2)     % this is for the marker 3 !
        hold on
       % plot(XX3,YY3,'.','MarkerSize',2)     % this is for the marker 3 !
        hold on
       % plot(XX4,YY4,'.','MarkerSize',2)     % this is for the marker 3 !
        %
        plot(XXX,YYY,'.','MarkerSize',2)    % this is for the inclusion !
        xlim([-3 3])
        ylim([-3 3])
        pbaspect([1 1 1])
        F(i)= getframe(f);
        %clf;
      
         parfor rrr=1:N                                % marker 1 bringing from general coordinate system to new rotated ellipsoid system
         newpos_plot(:,rrr)= q*newpos_plot(:,rrr);
         end
        %{
         parfor rrr=1:N1                                % marker 2 bringing from general coordinate system to new rotated ellipsoid system
         newpos_plot1(:,rrr)= q*newpos_plot1(:,rrr);
         end
         
         parfor rrr=1:N2                                % marker 3 bringing from general coordinate system to new rotated ellipsoid system
         newpos_plot2(:,rrr)= q*newpos_plot2(:,rrr);
         end
         
         parfor rrr=1:N3                                % marker 3 bringing from general coordinate system to new rotated ellipsoid system
         newpos_plot3(:,rrr)= q*newpos_plot3(:,rrr);
         end
         parfor rrr=1:N4                                % marker 3 bringing from general coordinate system to new rotated ellipsoid system
          newpos_plot4(:,rrr)= q*newpos_plot4(:,rrr);
         end
         %}
             
        %{     
        parfor rrr=1:NN                        % inclusion --bringing from general coordinate system to new rotated ellipsoid system
           epp(:,rrr)= q*epp1(:,rrr);
        end 
        %}
        %updating for next step  
   
         ep= newpos_plot;
         %{
         ep_1= newpos_plot1;
         ep_2= newpos_plot2;
         ep_3= newpos_plot3;
         ep_4=newpos_plot4;
         %}
     disp('Step completed');
     disp(i);
    toc 
end 
    

    
filename=sprintf('test2_FlankStruct_r=%d_ia%d_AR=%1.1f-%1.1f-%1.1f_Wk=%1.1f_grid3-3_5_markers_fine.avi',r,rad2deg(pp),A(1),A(2),A(3),Wk);
v= VideoWriter(filename);
v.FrameRate =2;
open(v)
writeVideo(v,F)
close(v)

      
 

      
        