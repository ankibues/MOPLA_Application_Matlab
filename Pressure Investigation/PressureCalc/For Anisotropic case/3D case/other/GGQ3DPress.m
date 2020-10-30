function T= GGQ3DPress(f, sub, x, a, b, c, d, Alp, Bet, ww,Cm)
% Global Gaussian Quadrature

% Input: f is a string,'fgin', 'pgin' or 'fgex','pgex';
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%        x, if f is 'in', x is a 3*1 matrix, 3 semi-axes of an inclusion;
%           if f is 'ex', x is a 3*2 matrix, 1st column is 3 semi-axes of 
%                  an inclusion, 2nd column is an coordinates of an point;
%        a,b is the integrating range of theta; c,d is the integrating range of phi;
%        gp is the number of gaussian points;

% Output: T is a scalar

  
  ft = funIn(f,sub,x,Alp,Bet,Cm);
    
  ww  = reshape(ww,1,[]);
  k   = ww.*ft;
  T = 0.25*(d-c)*(b-a)*sum(k);

end


function fr = funIn(ff,sub,a,alpha,beta,Cm)

   theta = reshape(alpha,1,[]);
   phi   = reshape(beta,1,[]);
   n = numel(theta);
  
   xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];
   aaa = a(1)*a(2)*a(3)/(4*pi);
   
   CC = Avinv(Cm,xi);
   
   aa = repmat(a,1,n);
   xia  = aa.*xi;
   xiasq = xia.^2;
   nxia = sqrt(sum(xiasq));
   rho = aaa*sin(phi)./nxia.^3;
   
   if strcmp(ff,'fgin')
       del = Delta(sub(1),sub(3));
       delta = repmat(del,1,n);
       fr = rho.*xi(sub(4),:).*xi(sub(2),:).*(delta - xi(sub(1),:).*xi(sub(3),:)); 
   elseif strcmp(ff,'pgin')
       fr = -rho.*CC(sub(1),:).*xi(sub(2),:);
   else
       fprintf('need to choose a function (fgin or pgin)');
   end
   
end

function re = Delta(m,n)
  % kroneckerDelta
  if m == n
      re = 1;
  else
      re = 0;
  end
end

function CC = Avinv(Cm,xi)
  [~,N]=size(xi);
  
Aij= zeros(3,3,N);



   for j=1:3
       for i=1:3
             %%%%%%% discuss this part
          Aij(i,j,:)= Cm(i,1,j,1)*xi(1,:).*xi(1,:) + Cm(i,1,j,2)*xi(1,:).*xi(2,:) + Cm(i,1,j,3)*xi(1,:).*xi(3,:)...
                        + Cm(i,2,j,1)*xi(2,:).*xi(1,:) + Cm(i,2,j,2)*xi(2,:).*xi(2,:) + Cm(i,2,j,3,:)*xi(2,:).*xi(3,:)...
                        + Cm(i,3,j,1)*(xi(3,:).*xi(1,:)) + Cm(i,3,j,2)*xi(3,:).*xi(2,:) + Cm(i,3,j,3)*xi(3,:).*xi(3,:);
       end
   end
      
 %%%%%--------- defining the Av matrix.
   a11 = Aij(1,1,:);
   a12 = Aij(1,2,:);
   a13 = Aij(1,3,:);
   Av(1,4,:)= xi(1,:); a14= Av(1,4,:);
   
   a21 = Aij(2,1,:);
   a22 = Aij(2,2,:);
   a23 = Aij(2,3,:);
   Av(2,4,:)= xi(2,:); a24= Av(2,4,:); 
   a31 = Aij(3,1,:);
   a32 = Aij(3,2,:);
   a33 = Aij(3,3,:);
   Av(3,4,:)= xi(3,:); a34= Av(3,4,:);
   Av(4,1,:,:)= xi(1,:);
   Av(4,2,:,:)= xi(2,:);
   Av(4,3,:,:)= xi(3,:);
   Av(4,4,:,:)= 0;
   a41= Av(4,1,:);     a42= Av(4,2,:);     a43= Av(4,3,:);    a44= Av(4,4,:);
 %%%%%-----following is a quick way to calculate inverse of the Av matrix

detA = a11.*a22.*a33.*a44 + a11.*a23.*a34.*a42 + a11.*a24.*a32.*a43...
      +a12.*a21.*a34.*a43 + a12.*a23.*a31.*a44 + a12.*a24.*a33.*a41...
      +a13.*a21.*a32.*a44 + a13.*a22.*a34.*a41 + a13.*a24.*a31.*a42...
      +a14.*a21.*a33.*a42 + a14.*a22.*a31.*a43 + a14.*a23.*a32.*a41...
      -a11.*a22.*a34.*a43 - a11.*a23.*a32.*a44 - a11.*a24.*a33.*a42...
      -a12.*a21.*a33.*a44 - a12.*a23.*a34.*a41 - a12.*a24.*a31.*a43...
      -a13.*a21.*a34.*a42 - a13.*a22.*a31.*a44 - a13.*a24.*a32.*a41...
      -a14.*a21.*a32.*a43 - a14.*a22.*a33.*a41 - a14.*a23.*a31.*a42 ;

  
A(1,1,:) = (a22.*a33.*a44 + a23.*a34.*a42 + a24.*a32.*a43 -a22.*a34.*a43-a23.*a32.*a44- a24.*a33.*a42)./detA ;
A(1,2,:) = (a12.*a34.*a43 + a13.*a32.*a44 + a14.*a33.*a42 -a12.*a33.*a44-a13.*a34.*a42- a14.*a32.*a43)./detA ;
A(1,3,:) = (a12.*a23.*a44 + a13.*a24.*a42 + a14.*a22.*a43 -a12.*a24.*a43-a13.*a22.*a44- a14.*a23.*a42)./detA ;
A(1,4,:) = (a12.*a24.*a33 + a13.*a22.*a34 + a14.*a23.*a32 -a12.*a23.*a34-a13.*a24.*a32- a14.*a22.*a33)./detA ;
A(2,1,:) = (a21.*a34.*a43 + a23.*a31.*a44 + a24.*a33.*a41 -a21.*a33.*a44-a23.*a34.*a41- a24.*a31.*a43)./detA ;
A(2,2,:) = (a11.*a33.*a44 + a13.*a34.*a41 + a14.*a31.*a43 -a11.*a34.*a43-a13.*a31.*a44- a14.*a33.*a41)./detA ;
A(2,3,:) = (a11.*a24.*a43 + a13.*a21.*a44 + a14.*a23.*a41 -a11.*a23.*a44-a13.*a24.*a41- a14.*a21.*a43)./detA ;
A(2,4,:) = (a11.*a23.*a34 + a13.*a24.*a31 + a14.*a21.*a33 -a11.*a24.*a33-a13.*a21.*a34- a14.*a23.*a31)./detA ;
A(3,1,:) = (a21.*a32.*a44 + a22.*a34.*a41 + a24.*a31.*a42 -a21.*a34.*a42-a22.*a31.*a44- a24.*a32.*a41)./detA ;
A(3,2,:) = (a11.*a34.*a42 + a12.*a31.*a44 + a14.*a32.*a41 -a11.*a32.*a44-a12.*a34.*a41- a14.*a31.*a42)./detA ;
A(3,3,:) = (a11.*a22.*a44 + a12.*a24.*a41 + a14.*a21.*a42 -a11.*a24.*a42-a12.*a21.*a44- a14.*a22.*a41)./detA ;
A(3,4,:) = (a11.*a24.*a32 + a12.*a21.*a34 + a14.*a22.*a31 -a11.*a22.*a34-a12.*a24.*a31- a14.*a21.*a32)./detA ;
A(4,1,:) = (a21.*a33.*a42 + a22.*a31.*a43 + a23.*a32.*a41 -a21.*a32.*a43-a22.*a33.*a41- a23.*a31.*a42)./detA ;
A(4,2,:) = (a11.*a32.*a43 + a12.*a33.*a41 + a13.*a31.*a42 -a11.*a33.*a42-a12.*a31.*a43- a13.*a32.*a41)./detA ;
A(4,3,:) = (a11.*a23.*a42 + a12.*a21.*a43 + a13.*a22.*a41 -a11.*a22.*a43-a12.*a23.*a41- a13.*a21.*a42)./detA ;
A(4,4,:) = (a11.*a22.*a33 + a12.*a23.*a31 + a13.*a21.*a32 -a11.*a23.*a32-a12.*a21.*a33- a13.*a22.*a31)./detA ;
  
% A is (Av)^-1 as in Eq. B4 Jiang (2014)

  

CC= [A(4,1,:);A(4,2,:);A(4,3,:)];
 
 
  end