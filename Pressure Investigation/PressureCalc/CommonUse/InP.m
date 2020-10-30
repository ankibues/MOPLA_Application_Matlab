function [h,u2,inp] = InP(x,D,Nm,Nc,Jd,Js,b)
 % This function is to get pressure inside the clast and other values needed in the
 % calculation of pressure outside the clast
 
  q = x(1:3,1:3);
  d = q*D*q';
  a = x(:,4);
  
  [s,p] = SnpIn(a,Jd,Js);
  invS    = Inverse(s,b);
  h       = Contract(2*Jd,invS);
  p_in    = R_Multiply(p,h);
  E       = Ed(Nm,Nc,x(1,5),s,d,x(2,5),b,Jd);
  dE      = E(1:3,1:3);
  u2      = dE-d;
  inp     = contract1(p_in,u2);
end