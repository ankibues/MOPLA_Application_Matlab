% HSZM_Cylinder_isotropic_newtonian

% Summary: The analytic solutions of the Eshelby tnsors for both interior
% and exterior fields are used.
 % Jan-14-2019
 % clear all variables, Comment Windows and figures
   clear;
   clc;
 % clf;
 % Input parameters:
 % imposed far-field flow L (exterior coordinate)
 % xy plane-straining general shearing
 L = [.5   1   0;...
      0   -.5   0;...
      0   0   0];
 % decompose the far-field flow L into a strain rate tensor D and a vorticity 
 % tensor W, Eqn(3) in Jiang(2007a)     
 D     = 0.5 * (L + L');
 W     = 0.5 * (L - L');  
 %epsilonII = 0.5;
 % Cylinder inclusion
 a = [100; 20; 1];
 % the viscosity ratio at matrix strain rate, eta_n/eta_matrix
 r = 0.01;
 % orientation
 theta = deg2rad(0);
 q = [0 0 1;cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0]; 
 %q = [0 0 1;1 0 0;0 1 0];
%  % stress exponent of the matrix
%   Nm    = 1; 
% % stress exponent of the ellipsoid 
%   Ne    = 1;  
  
 % time increment
 tincr = 0.1;
 % total steps
 step = 5;
 
 % generate 4th-order identity tensors   
 [Jd, ~, Ja, ~] = FourIdentity();   
 % stiffness tensor of the matrix, assuming it's isotropic here, in interior
 % coordinate.
 Cm    = 2*Jd;
 % stiffness tensor of the ellipsoid, assuming it's a planar anisotropic
 % material with the anisotropy plane normal to a3-axis,in interior coordinate. 
 Ce    = 2*r*Jd;  
 
 % position coordinate in global coordinate        
  xgv     = -25:1:25;           % grid vector: x'axis,a1
  ygv     = -5:1:5;                   % grid vector: y'axis,a2
  zgv     = 0;             % grid vector: z'axis,a3
  [X,Y,Z] = meshgrid(xgv,ygv,zgv);
  pc      = [reshape(X,1,[]);reshape(Y,1,[]);reshape(Z,1,[])]; 
  
  num     = numel(X);
% bring the mesh to clast's coordinate
  ppc    = pc;
  for rrr=1:num
      ppc(:,rrr) = q * pc(:,rrr);
  end
  
% Check if the points are inside the ellipsoid or not
%   ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
%   x_ex  = X(ind);
%   y_ex  = Y(ind);
%   z_ex  = Z(ind);
%   ep    = cat(1,x_ex',y_ex',z_ex');
   ind   = (ppc(1,:)./a(1)).^2 + (ppc(2,:)./a(2)).^2 + (ppc(3,:)./a(3)).^2 > 1;
   ep    = ppc(:,ind);
   [~,n] = size(ep);  
   
   
   %% Interior Fields 
%--------------------------------------------------------------------------
% Eshelby Tensors (S,PI) for Interior points
% Elastic Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el,PI_el] = SnP(a);

% Viscous Eshelby Tensors (S,PI,LAMDA) for Interior points
% LAMBDA
  LAMBDA   = zeros(3,3);
  for i=1:3
      LAMBDA(i,i) = -1/3* (S_el(i,i,1,1)+ S_el(i,i,2,2)+ S_el(i,i,3,3)); 
  end
% S  
  S       = S_el;
  for i=1:3
      for j=1:3
          S(i,i,j,j) = LAMBDA(i,i)+ S_el(i,i,j,j);
      end
  end
% PI  
  PI      = PI_el;
  
% Interior fields (e,w,s) in the clast's coordinate system
% describe D,W in the clast's coordinate system 
  D_bar   = q * D * q';
  W_bar   = q * W * q';
% e: interior field strain-rate
  invS    = FourTensorInv(S);
  H       = Contract(Cm,(invS - Jd));   %Hill's constraint tensor
  t1      = FourTensorInv(H + Ce);
  A       = Contract(t1,(H + Cm));      %strain-rate partitioning tensor
  e       = Multiply(A, D_bar); 
  
  %[e1,C_clst] = Ed(Nm,Ne,Cm,Ce,r,S,D_bar,epsilonII,Jd,q);
  
% s: interior field stress
  s        = Multiply(Ce, e);
  s_m      = q' * s * q;
  [VV,DD]  = eig(s_m);
  [~,indd]= sort(diag(DD));
  svec     = VV(:,indd);
 
% w: interior field vorticity
  invS     = FourTensorInv(S);
  de       = e - D_bar;                  
  t4       = Multiply(Contract(PI, invS), de);
  w        = t4 + W_bar;

% l: interior flow in ellipsoid coordinate
  l        = e + w;  
% f: interior position gradient tensor
  f        = expm(l*tincr*step);
  %transfer f back to global coordinate
  f_m      = q' * f * q;
% the eigenvalues and aigenvectors of finite strain ellipsoid
% eval=[e1;e2;e3], e1>e2>e3; evec=[v1 v2 v3] column vectors
  [eval,evec]=finite_strain_ellipsoid(f_m);
  fvec1 = evec(:,1);
  fvec3 = evec(:,3);
  
% bulk finite strain ellipsoid
  FF        = expm(L*tincr*step);
  [Feval,evec]=finite_strain_ellipsoid(FF);
  Fvec1 = evec(:,1);
  Fvec3 = evec(:,3);
  Sigma = Multiply(Cm,D);
  [VV,DD]  = eig(Sigma);
  [~,indd]= sort(diag(DD));
  Svec     = VV(:,indd);
 %% Exterior Fields  
%--------------------------------------------------------------------------  
  A = [a(2);a(3)];
  
  F = repmat(eye(3),1,1,n);
  F_evl  = repmat(F,1,1,1,step);
  ep_evl = zeros(3,n,step);
  vel_f  = zeros(2,n);
%   A1       = zeros(3,n,step);
%   e1       = zeros(3,n,step);
%   e3       = zeros(3,n,step);
%   ang1     = zeros(2,n,step);
%   ang3     = zeros(2,n,step);
%   r1       = zeros(n,step);
%   r3       = zeros(n,step);

  for tt=1:step
      tic
      ind1 = (ep(1,:)./a(1)).^2 + (ep(2,:)./a(2)).^2 +...
             (ep(3,:)./a(3)).^2 > 1;
      for ii=1:n
          if ind1(ii)==1
              %points outside the inclusion
              %bring the points to Elliptic coordinate
              [xi,eta]= CartesiantoElliptic(ep(2,ii),ep(3,ii),A);
              % choosing the known points for velocity field
              % these known points in elliptic coordinates
              [xi0,~]= CartesiantoElliptic(A(1),0,A);    
              eta0= real(eta);
              % these known points in cartesian coordinates, 
              % so that velocity field can be calculated.
              [x0,y0]= ElliptictoCartesian(xi0,eta0,A);
              % reference velocity field
              v0 = l*[0;x0;y0];
              
              % the interval of integration 
              kk  = .001;
              XI  = xi0:kk:xi;
              ETA = eta0*ones(size(XI));
              % bring them to cartesian coordinate system
              [x,y]= ElliptictoCartesian(XI,ETA,A);

              grid = [zeros(size(x)); x; y];
              [~,N1]= size(grid);
              % The exterior Eshelby tensors(S_Ex PI_Ex) 
              [S_Ex,PI_Ex] = SnPI_Ex(a,grid);

              % Exterior fields 
              e_Ex = zeros(3,3,N1);
              w_Ex = zeros(3,3,N1);
              l_Ex = zeros(3,3,N1);
              
              parfor rr=1:N1
                  
                  % e_Ex: exterior field strain-rate in clast's coordinate
                  v1           = Contract(S_Ex(:,:,:,:,rr), invS);
                  v2           = Multiply(v1, de);
                  e_Ex(:,:,rr) = v2 + D_bar;
                  % w_Ex: exterior field vorticity in clast's coordinate
                  v3           = Contract(PI_Ex(:,:,:,:,rr), invS);
                  v4           = Multiply(v3, de);
                  w_Ex(:,:,rr) = v4 + W_bar; 
                  % l_Ex: exterior field flow in clast's coordinate
                  l_Ex(:,:,rr)  = e_Ex(:,:,rr) + w_Ex(:,:,rr);
              
              end
              % calculate the intervals.
              % Note: the invertals are in elliptic system but not equal in carterian
              % system. So this process has to be done.
              grid1 = [zeros(3,1),grid(:,1:N1-1)];
              delta = grid - grid1;
              delta(:,1)=[];
              deltaX = delta(2,:);
              deltaY = delta(3,:);
            
              % update the position gradiant tensor.
              F(:,:,ii) = (eye(3) + l_Ex(:,:,N1)*tincr)* F(:,:,ii);
            
              % update the velocity field.
              l_Ex(:,:,N1) = [];
              V2= v0(2) + sum(deltaX.*(squeeze(l_Ex(2,2,:)))',2) + sum(deltaY.*(squeeze(l_Ex(2,3,:)))',2);
              V3= v0(3) + sum(deltaX.*(squeeze(l_Ex(3,2,:)))',2) + sum(deltaY.*(squeeze(l_Ex(3,3,:)))',2);
              
              vel_f(:,ii) = [V2;V3];
              
              % update the coordinate.
              ep(2,ii) = ep(2,ii)+ V2*tincr;
              ep(3,ii) = ep(3,ii)+ V3*tincr;
          else
              % update the position gradiant tensor.
              F(:,:,ii) = (eye(3) + l*tincr)* F(:,:,ii);
              % update the velocity field.
              vv = l*ep(:,ii);
              % update the coordinate.
              ep(2,ii) = ep(2,ii)+ vv(2)*tincr;
              ep(3,ii) = ep(3,ii)+ vv(3)*tincr;
              
              vel_f(:,ii) = [vv(2);vv(3)];
          end
          
          % Transfer point and the position grediant tensor F into 
          % global coordinate 
          ep_evl(:,ii,tt)  = q' * ep(:,ii);
          F_evl(:,:,ii,tt) = q' * F(:,:,ii) * q;
         
          
      end
      
      XX=['Step ',num2str(tt),' completed'];
      disp(XX);
      toc
  end
%  

%  
F_final = F_evl(:,:,:,end);
Fe1     = zeros(3,n);
Fe3     = zeros(3,n);
for rr=1:n
   [eval,evec]=finite_strain_ellipsoid(F_final(:,:,rr));
   Fe1(:,rr) = evec(:,1);
   Fe3(:,rr) = evec(:,3);
end
% colorVec = jet(n);
%  syms x y
%   circle1 = x^2/25 + y^2 == 1;
%   fimplicit(circle1)  
%   axis equal 
%   hold on
%   for ii=1:n
%       %plot(pc_evl(1,ii,1),pc_evl(2,ii,1),'o','Color',colorVec(ii,:))
%       %plot(squeeze(pc_evl(1,ii,2:end-1)),squeeze(pc_evl(2,ii,2:end-1)),'.','Color',colorVec(ii,:))
%       %plot(ep_evl(1,ii,end),ep_evl(2,ii,end),'*','Color',colorVec(ii,:))
%       quiver(0,0,fvec1(1),fvec1(2),'m')
%       quiver(0,6,Fvec1(1),Fvec1(2),'r')
%       quiver(0,0,svec(1,1),svec(2,1),'b')
%       quiver(0,6,Svec(1,1),Svec(2,1),'g')
%       quiver(ep_evl(1,ii,end),ep_evl(2,ii,end),Fe1(1,ii),Fe1(2,ii),'m')
%      
%   end
%  hold off 
XX = squeeze(ep_evl(1,:,end));
YY = squeeze(ep_evl(2,:,end));
FX = Fe1(1,:);
FY = Fe1(2,:);

figure('Name','the long axes of the finite strain ellipsoids'); 
 syms x y
 circle1 = (x/a(2))^2 + y^2 == 1;
 fimplicit(circle1) 
 hold on
q=quiver(XX, YY, FX, FY,'b');
q.ShowArrowHead = 'off';
q.AutoScale = 'on';
q.AlignVertexCenters = 'on';
hold on
quiver(0,0,fvec1(1),fvec1(2),'m')
hold on
plot(squeeze(ep_evl(1,:,end)),squeeze(ep_evl(2,:,end)),'.')
xlabel('x')
ylabel('y')
axis equal