function r = GLeSH(f, sub, x, theta, phi, w)
% Global Lebedev Sphrical Harmonics

% Input: f is a string, 'fgin', 'pgin' or 'fgex','pgex';
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%        x, if f is 'in', x is a 3*1 matrix, 3 semi-axes of an inclusion;
%           if f is 'ex', x is a 3*2 matrix, 1st column is 3 semi-axes of 
%                  an inclusion, 2nd column is an coordinates of an point;
%        ld is the number of Lebedev points;
%        choice:6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,
%      770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810
  

% Output: r is a scalar

  
  if strfind(f,'in')
      ft = funIn(f,sub,x,theta,phi);
  elseif strfind(f,'ex')
      ft = funEx(f,sub,x(:,1),x(:,2),theta,phi);
  else
      fprintf('please choose between in and ex');
  end
 
  w   = reshape(w,1,[]);
  k   = w.*ft;
  r  = 4*pi*sum(k);

end

function fr = funIn(ff,sub,a,alpha,beta)

   theta = reshape(alpha,1,[]);
   phi   = reshape(beta,1,[]);
   n = numel(theta);
  
   xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];
   aaa = a(1)*a(2)*a(3)/(4*pi);
   
   aa = repmat(a,1,n);
   xia  = aa.*xi;
   xiasq = xia.^2;
   nxia = sqrt(sum(xiasq));
   rho = aaa./nxia.^3;
   
   if strcmp(ff,'fgin')
       del = Delta(sub(1),sub(3));
       delta = repmat(del,1,n);
       fr = rho.*xi(sub(4),:).*xi(sub(2),:).*(delta - xi(sub(1),:).*xi(sub(3),:));
   elseif strcmp(ff,'pgin')
       fr = -rho.*xi(sub(1),:).*xi(sub(2),:);
   else
       fprintf('need to choose a function (fgin or pgin)');
   end
   
end

function fr = funEx(ff,sub,a,x,Alp,Bet)
  
  theta = reshape(Alp,1,[]);
  phi   = reshape(Bet,1,[]);
  n = numel(theta);
  
  xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];
  aaa     = a(1)*a(2)*a(3)/(4*pi);
  
  aa = repmat(a,1,n);
  xx = repmat(x,1,n);
  
  epsa  = aa.*xi;
  A     = xi./aa;
  r     = xx-epsa;
  
  rsq = r.^2;
  nr = sqrt(sum(rsq));
  
  temp = aaa./nr.^3;
  rho(1,:) = temp;
  rho(2,:) = temp;
  rho(3,:) = temp;
  rhoA = rho.*A;
  
  if strcmp(ff,'fgex')
      fr = 0.5*rhoA(sub(2),:).*(-Delta(sub(1),sub(3))*r(sub(4),:)...
               +Delta(sub(1),sub(4))*r(sub(3),:)+Delta(sub(3),sub(4))*r(sub(1),:)-...
               3*(r(sub(1),:).*r(sub(3),:).*r(sub(4),:))./(nr.^2));
  elseif strcmp(ff,'pgex')  
      fr = r(sub(1),:).*rhoA(sub(2),:);
  else
      fprintf('need to choose a function (fgex or pgex)');
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