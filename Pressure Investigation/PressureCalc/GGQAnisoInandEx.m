function T = GGQAnisoInandEx(f, sub, x, Alp, ww,phi,m,Nn)
% Global Gaussian Quadrature

% Input: f is a string,'fgin', 'pgin' or 'fgex','pgex';
%        sub is a 1*4 matrix, 4 subscripts of the 4th order tensor T;
%        x, if f is 'in', x is a 2*1 matrix, 2 semi-axes of an inclusion;
%           if f is 'ex', x is a 2*2 matrix, 1st column is 2 semi-axes of 
%                  an inclusion, 2nd column is an coordinates of an point;
%        a,b is the integrating range of psi
%        gp is the number of gaussian points;

% Output: T is a scalar

  
  if strfind(f,'in')
      A= x(1);
      B= x(2);
      ft = funIn(f,sub,x,Alp,phi,m,Nn);
  elseif strfind(f,'ex')
      A= x(1,1);
      B= x(2,1);
      ft = funEx(f,sub,x(:,1),x(:,2),Alp,phi,m,Nn);
  else
      fprintf('please choose between in and ex');
  end
  
  
  k   = ww'.*ft;
  T =  A*B*sum(k)/pi;

end


function fr = funIn(ff,sub,a,alpha,phi,m,Nn)
  
   psi = alpha';
   
   n = numel(psi);
   z = [cos(psi); sin(psi)];
   epsilon = (m*(sin(2*(phi + psi))).^2 + (cos(2*(phi + psi))).^2).^(-1);
   epsilon_rep = repmat(epsilon, 2,1)/2;
   C = epsilon_rep.*[(m+1)*cos(psi)-(m-1)*cos(4*phi + 3*psi); (m+1)*sin(psi)+(m-1)*sin(4*phi + 3*psi)];
           
   B = 1./((a(1)*cos(psi)).^2 + (a(2)*sin(psi)).^2);
   
   if strcmp(ff,'fgin')
       del = Delta(sub(1),sub(3));
       delta = repmat(del,1,n);
       fr = B.*z(sub(2),:).*z(sub(4),:).*(delta - (z(sub(1),:).*z(sub(3),:)))*(m/Nn).*epsilon; 
   elseif strcmp(ff,'pgin')
       fr = -B.*C(sub(1),:).*z(sub(2),:);
   else
       fprintf('need to choose a function (fgin or pgin)');
   end
   
end

function fr = funEx(ff,sub,a,x,Alp,phi,m,Nn)
  
  psi = Alp';
  n = numel(psi);
  
  z = [cos(psi); sin(psi)];
  aaa = (m/Nn);
  epsilon = (m*(sin(2*(phi + psi))).^2 + (cos(2*(phi + psi))).^2).^(-1);
  epsilon_rep= repmat(epsilon, 2,1)/2;
  C = epsilon_rep.*[(m+1)*cos(psi)-(m-1)*cos(4*phi + 3*psi); (m+1)*sin(psi)+(m-1)*sin(4*phi + 3*psi)];
  
  B = 1./((a(1)*cos(psi)).^2 + (a(2)*sin(psi)).^2);
  k= zeros(1,numel(psi));
  for o=1:numel(psi)
      k(1,o) = real(1 - ((1i*(abs(dot(x,z(:,o)))))/(((a(1)*z(1,o))^2 + (a(2)*z(2,o))^2 - (dot(x,z(:,o)))^2  )^0.5)));
  end
      
    
  if strcmp(ff,'fgex')
      del = Delta(sub(1),sub(3));
      delta = repmat(del,1,n);
      fr = aaa*B.*k.*z(sub(2),:).*z(sub(4),:).*(delta - z(sub(1),:).*z(sub(3),:)).*epsilon;
  elseif strcmp(ff,'pgex')  
      fr = -B.*k.*C(sub(1),:).*z(sub(2),:);
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
