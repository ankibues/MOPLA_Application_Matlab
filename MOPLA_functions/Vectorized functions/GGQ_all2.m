function T = GGQ_all2(sub,x, a, b, c, d, alpha, beta, ww,Cm)
% Global Gaussian Quadrature solution for 'N' inclusions

% Input: f is a string,'fgin', 'pgin' or 'fgex','pgex';
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%          x is a 3*N matrix, 3 semi-axes of 'N' inclusions;
%        Cm is the stiffness tensor 3*3*3*3*N in the frame of reference of 'N' number of ellipsoids.           
%        a,b is the integrating range of theta; c,d is the integrating range of phi;
%        gp is the number of gaussian points;

% Output: T is a scalar

% ff should be N*10000 matrix.
   [~,N]=size(x);
   theta = reshape(alpha,1,[]);
   phi   = reshape(beta,1,[]);
   n = numel(theta);
  
   xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];
   aaa = x(1,:).*x(2,:).*x(3,:)./(4*pi);
   aaa=aaa';
   
   
   %%%%%%%%%%%%%%%%%%%%%%%--- Denominator term of the integral (eq. B4,
   %%%%%%%%%%%%%%%%%%%%%%%-----Jiang 2014)
   aa = repmat(x,1,1,n);
   A1= squeeze(aa(1,:,:));
   XIA1 = A1.*repmat(xi(1,:),N,1);
   A2= squeeze(aa(2,:,:));
   XIA2 = A2.*repmat(xi(2,:),N,1);
   A3= squeeze(aa(3,:,:));
   XIA3 = A3.*repmat(xi(3,:),N,1);
   
   denom= ((XIA1).^2 + (XIA2).^2 +(XIA3).^2).^(3/2); 
   %%%%%%%%%%%%%%%%%%%%%%%------Numerator term of the integral (eq. B4,
   %%%%%%%%%%%%%%%%%%%%%%%------Jiang 2014)

  numer= repmat(xi(sub(2),:).*xi(sub(4),:).*sin(phi),N,1).* Avinv(sub(1),sub(3),Cm,xi);
  ft = numer./denom; 
  
  ww  = reshape(ww,1,[]);
  www= repmat(ww,N,1);
  k   = www.*ft;
  T = .25*(b-a)*(d-c)* aaa.*sum(k,2);%   
   
end


  function Avinvvv = Avinv(I,J,Cm,xi)
  [~,N]=size(xi);
  [~,~,~,~,n]= size(Cm);
Aij= zeros(3,3,n,N);



   for j=1:3
       for i=1:3
             %%%%%%% discuss this part
          Aij(i,j,:,:)= squeeze(Cm(i,1,j,1,:))*(xi(1,:).*xi(1,:)) + squeeze(Cm(i,1,j,2,:))*(xi(1,:).*xi(2,:)) + squeeze(Cm(i,1,j,3,:))*(xi(1,:).*xi(3,:))...
                      + squeeze(Cm(i,2,j,1,:))*(xi(2,:).*xi(1,:)) + squeeze(Cm(i,2,j,2,:))*(xi(2,:).*xi(2,:)) + squeeze(Cm(i,2,j,3,:))*(xi(2,:).*xi(3,:))...
                      + squeeze(Cm(i,3,j,1,:))*(xi(3,:).*xi(1,:)) + squeeze(Cm(i,3,j,2,:))*(xi(3,:).*xi(2,:)) + squeeze(Cm(i,3,j,3,:))*(xi(3,:).*xi(3,:));
       end
   end
      
 %%%%%--------- defining the Av matrix.
   a11 = Aij(1,1,:,:);
   a12 = Aij(1,2,:,:);
   a13 = Aij(1,3,:,:);
   Av(1,4,:,:)= repmat(xi(1,:),n,1); a14= Av(1,4,:,:);
   
   a21 = Aij(2,1,:,:);
   a22 = Aij(2,2,:,:);
   a23 = Aij(2,3,:,:);
   Av(2,4,:,:)= repmat(xi(2,:),n,1); a24= Av(2,4,:,:); 
   a31 = Aij(3,1,:,:);
   a32 = Aij(3,2,:,:);
   a33 = Aij(3,3,:,:);
   Av(3,4,:,:)= repmat(xi(3,:),n,1); a34= Av(3,4,:,:);
   Av(4,1,:,:)= repmat(xi(1,:),n,1);
   Av(4,2,:,:)= repmat(xi(2,:),n,1);
   Av(4,3,:,:)= repmat(xi(3,:),n,1);
   Av(4,4,:,:)= 0;
   a41= Av(4,1,:,:);     a42= Av(4,2,:,:);     a43= Av(4,3,:,:);    a44= Av(4,4,:,:);
 %%%%%-----following is a quick way to calculate inverse of the Av matrix

detA = a11.*a22.*a33.*a44 + a11.*a23.*a34.*a42 + a11.*a24.*a32.*a43...
      +a12.*a21.*a34.*a43 + a12.*a23.*a31.*a44 + a12.*a24.*a33.*a41...
      +a13.*a21.*a32.*a44 + a13.*a22.*a34.*a41 + a13.*a24.*a31.*a42...
      +a14.*a21.*a33.*a42 + a14.*a22.*a31.*a43 + a14.*a23.*a32.*a41...
      -a11.*a22.*a34.*a43 - a11.*a23.*a32.*a44 - a11.*a24.*a33.*a42...
      -a12.*a21.*a33.*a44 - a12.*a23.*a34.*a41 - a12.*a24.*a31.*a43...
      -a13.*a21.*a34.*a42 - a13.*a22.*a31.*a44 - a13.*a24.*a32.*a41...
      -a14.*a21.*a32.*a43 - a14.*a22.*a33.*a41 - a14.*a23.*a31.*a42 ;

  
A(1,1,:,:) = (a22.*a33.*a44 + a23.*a34.*a42 + a24.*a32.*a43 -a22.*a34.*a43-a23.*a32.*a44- a24.*a33.*a42)./detA ;
A(1,2,:,:) = (a12.*a34.*a43 + a13.*a32.*a44 + a14.*a33.*a42 -a12.*a33.*a44-a13.*a34.*a42- a14.*a32.*a43)./detA ;
A(1,3,:,:) = (a12.*a23.*a44 + a13.*a24.*a42 + a14.*a22.*a43 -a12.*a24.*a43-a13.*a22.*a44- a14.*a23.*a42)./detA ;
A(1,4,:,:) = (a12.*a24.*a33 + a13.*a22.*a34 + a14.*a23.*a32 -a12.*a23.*a34-a13.*a24.*a32- a14.*a22.*a33)./detA ;
A(2,1,:,:) = (a21.*a34.*a43 + a23.*a31.*a44 + a24.*a33.*a41 -a21.*a33.*a44-a23.*a34.*a41- a24.*a31.*a43)./detA ;
A(2,2,:,:) = (a11.*a33.*a44 + a13.*a34.*a41 + a14.*a31.*a43 -a11.*a34.*a43-a13.*a31.*a44- a14.*a33.*a41)./detA ;
A(2,3,:,:) = (a11.*a24.*a43 + a13.*a21.*a44 + a14.*a23.*a41 -a11.*a23.*a44-a13.*a24.*a41- a14.*a21.*a43)./detA ;
A(2,4,:,:) = (a11.*a23.*a34 + a13.*a24.*a31 + a14.*a21.*a33 -a11.*a24.*a33-a13.*a21.*a34- a14.*a23.*a31)./detA ;
A(3,1,:,:) = (a21.*a32.*a44 + a22.*a34.*a41 + a24.*a31.*a42 -a21.*a34.*a42-a22.*a31.*a44- a24.*a32.*a41)./detA ;
A(3,2,:,:) = (a11.*a34.*a42 + a12.*a31.*a44 + a14.*a32.*a41 -a11.*a32.*a44-a12.*a34.*a41- a14.*a31.*a42)./detA ;
A(3,3,:,:) = (a11.*a22.*a44 + a12.*a24.*a41 + a14.*a21.*a42 -a11.*a24.*a42-a12.*a21.*a44- a14.*a22.*a41)./detA ;
A(3,4,:,:) = (a11.*a24.*a32 + a12.*a21.*a34 + a14.*a22.*a31 -a11.*a22.*a34-a12.*a24.*a31- a14.*a21.*a32)./detA ;
A(4,1,:,:) = (a21.*a33.*a42 + a22.*a31.*a43 + a23.*a32.*a41 -a21.*a32.*a43-a22.*a33.*a41- a23.*a31.*a42)./detA ;
A(4,2,:,:) = (a11.*a32.*a43 + a12.*a33.*a41 + a13.*a31.*a42 -a11.*a33.*a42-a12.*a31.*a43- a13.*a32.*a41)./detA ;
A(4,3,:,:) = (a11.*a23.*a42 + a12.*a21.*a43 + a13.*a22.*a41 -a11.*a22.*a43-a12.*a23.*a41- a13.*a21.*a42)./detA ;
A(4,4,:,:) = (a11.*a22.*a33 + a12.*a23.*a31 + a13.*a21.*a32 -a11.*a23.*a32-a12.*a21.*a33- a13.*a22.*a31)./detA ;
  
%Avinvv= Ainv(Av,n,N);

  

Avinvvv= reshape(A(I,J,:,:),n,N);
 
 
  end
%{
function A= Ainv(Av,n,N)
% calculate inverse of all Av at once without loop

A= zeros(4,4,n,N);








end
%}


