%function Single_marker_FS_3D(Wk,r,a,pp)
% Deformation of marker with exterior points around inclusion- modelling a
% flanking structure in 3D 
% Here, marker is a plane.

k=pwd;
 K= extractBefore(k,"Flanking Structure Investigation");
 addpath(genpath(K));
%parpool(32);
tic

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
   steps =1000;
   for ii=1:steps
       [rho,~,~]=Calc_Gamma(L,tincr,ii);
       if rho>2
           break
       end
   end
   steps=ii;
       
%  Power law coefficients for matrix and inclusion
   Nm=1;
   Nc=1; 
   r= .1;
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
   a= [1;.5;.1];
   A=a;% initial length of the ellipsoid axis
   pp= 3*pi/4;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
 
% Here, I define the markers in the external coordinate system

% 3D meshgrid in external coordinate
        xgv     = -1.5:.1:1.5;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = zeros(size(xgv));    % grid vector: y'axis,a2 ----since the points are in X-Z plane
        zgv     = -1.5:.1:1.5;  %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        
        [Xx,Yy,Zz]= meshgrid(xgv,ygv,zgv);
        [~,~,nnn]= size(Xx);
        k1= Xx(:);
        k2= Yy(:);
        k3= Zz(:);
        
        ep    = [k1';k2';k3'];

        [~,N]= size(ep);
        

        
   MarkerPoints= zeros(3,N,steps);     
    QA_steps= zeros(3,4,steps); 
      % bringing markers into interior coordinate system
  
    for rrr=1:N                               
        ep(:,rrr)= q*ep(:,rrr);
    end
    
   

    
% frames to capture in every step:
  
 % F(steps)= struct('cdata',[],'colormap',[]);
 
 %

for i=1:steps
   
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

     


    % Now, once all the exterior points are available , we calculate
    % velocity field and displacement field for them.
   %------------------------------------------------------ 
    % for first marker !
    vel_f= zeros(3,N);
    xxx= ep(1,:);
    yyy= ep(2,:);
    zzz= ep(3,:);
    
    % [~,N]=size(ep);
    %
    a1=a(1);
    a2=a(2);
    a3=a(3);
    
    parfor ir=1:N
        if (xxx(ir)/a1)^2 + (yyy(ir)/a2)^2 + (zzz(ir)/a3)^2 <= 1
            vel_f(:,ir)=ll*[xxx(ir);yyy(ir);zzz(ir)];
        else
            vel_f(:,ir)= Velocity_field_calc_general3D(xxx(ir),yyy(ir),zzz(ir),a,invS,Jd,d,u2,w,ll);
        end
        
    end
    
    disp_f= vel_f*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form
    
    kkk= zeros(1,N);
    newpos_plot= zeros(size(disp_f));
    newpos= ep + disp_f;  % this is the new position calculated for the marker points.

    
     
      
     for rrr=1:N                                % bringing to general coordinate system for plotting
         newpos_plot(:,rrr)= q'*newpos(:,rrr);
     end
         
   
     MarkerPoints(:,:,i)= newpos_plot;
     QA_steps(:,:,i)=[q a];
     
     %------------------------------------------------------------
 
    
        %{   
        % 3D meshgrid in internal coordinate -----interior points of
        % inclusion-------------------------------------------------------------------------
        [X1,Y1,Z1]= ellipsoid(0,0,0,a(1),a(2),a(3));
        XX1=X1(:);
        YY1=Y1(:);
        ZZ1=Z1(:);
        pts= [XX1;YY1;ZZ1];    % these are the points on the surface of ellipsoid in ellipsoids coordinate system
        [~,NN]= size(pts);
        
        for rrr=1:NN                                % bringing pts to general coordinate system
            pts(:,rrr)= q'*pts(:,rrr);
        end 
        
        XX1= pts(1,:);
        YY1= pts(2,:);
        ZZ1= pts(3,:);
        
        % getting points in the form of meshgrid
        X_mesh= reshape(XX1,[21,21]);
        Y_mesh= reshape(YY1,[21,21]);
        Z_mesh= reshape(ZZ1,[21,21]);
        %}
        % this can be plotted using the surf command
        if i==steps
          save('MarkerPos3D_finegrid.mat','MarkerPoints','QA_steps','newpos_plot','a','q','steps','Wk','r'); 
        end
        
        qanew= QL(a,q,theta,dE,tincr);  
         a=qanew(:,4);
         q=qanew(1:3,1:3);
        
       % plotting part for each step
       % bringing points back to meshgrid for plotting
       
      % XX= reshape(newpos_plot(1,:),nnn,nnn,nnn);
       %YY= reshape(newpos_plot(2,:),nnn,nnn,nnn);
     %  ZZ= reshape(newpos_plot(3,:),nnn,nnn,nnn);
       
       
       
   
       
      
      % inclusion
       % f= figure('visible','on');
      
       % surf(XX,YY,ZZ)      % this is for the marker 1 !
      %  hold on      
       % surf(X_mesh,Y_mesh,Z_mesh)    % this is for the inclusion !
       % xlim([-3 3])
       % ylim([-3 3])
        %zlim([-3 3])
       % pbaspect([1 1 1])
       %{ 
        F(i)= getframe(f);
        [rho,~,~]=Calc_Gamma(L,tincr,i);
        if rho>2.4 && rho<2.5
            filename1=sprintf('FlankStruct_r=%d_ia%d_AR=%1.1f-%1.1f-%1.1f_Wk=%1.1f_rho=%1.1f.jpg',r,rad2deg(pp),A(1),A(2),A(3),Wk,rho);
            saveas(f,filename1);
        end
        if rho>5
            filename1=sprintf('FlankStruct_r=%d_ia%d_AR=%1.1f-%1.1f-%1.1f_Wk=%1.1f_rho=%1.1f.jpg',r,rad2deg(pp),A(1),A(2),A(3),Wk,rho);
            saveas(f,filename1);
            save('MarkerPos1.mat','XX','YY','a','q');
            break
        end
     %}              
        %clf;
         
         for rrr=1:N                                % marker 1 bringing from general coordinate system to new rotated ellipsoid system
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
     
     disp('rho=');
     disp(rho)
     disp('Step completed');
     disp(i);
    
    
end 
    

    toc
%{    
filename=sprintf('22FlankStruct_r=%d_ia%d_AR=%1.1f-%1.1f-%1.1f_Wk=%1.1f_rho=%1.1f.avi',r,rad2deg(pp),A(1),A(2),A(3),Wk,rho);
v= VideoWriter(filename);
v.FrameRate =2;
open(v)
writeVideo(v,F)
close(v)
%end
%}
      
      
        