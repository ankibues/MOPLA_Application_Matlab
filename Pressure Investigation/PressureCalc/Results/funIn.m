function T = funIn(ff,sub,a,alpha,phi,m,Nn,ww)
  
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
   
   k   = ww'.*fr;
  T =  a(1)*a(2)*sum(k)/pi;
end