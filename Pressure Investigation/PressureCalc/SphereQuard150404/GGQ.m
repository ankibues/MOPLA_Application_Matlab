function T = GGQ(f, sub, x, a, b, c, d, Alp, Bet, ww)
% Global Gaussian Quadrature

% Input: f is a string,'fgin', 'pgin' or 'fgex','pgex';
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%        x, if f is 'in', x is a 3*1 matrix, 3 semi-axes of an inclusion;
%           if f is 'ex', x is a 3*2 matrix, 1st column is 3 semi-axes of 
%                  an inclusion, 2nd column is an coordinates of an point;
%        a,b is the integrating range of theta; c,d is the integrating range of phi;
%        gp is the number of gaussian points;

% Output: T is a scalar

  
  if strfind(f,'in')
      ft = funIn(f,sub,x,Alp,Bet);
  elseif strfind(f,'ex')
      ft = funEx(f,sub,x(:,1),x(:,2),Alp,Bet);
  else
      fprintf('please choose between in and ex');
  end
  
  ww  = reshape(ww,1,[]);
  k   = ww.*ft;
  T = 0.25*(d-c)*(b-a)*sum(k);

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
   rho = aaa*sin(phi)./nxia.^3;
   
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
  
  temp = aaa*sin(phi)./nr.^3;
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
