function result = fdE(n,rr,S,b,Jd)
  
  a      = Inverse(Jd+((n*rr)-1)*S,b);    % eq. 16 of Jiang 2012b(deformable ellipsoid in power law viscous materials)
  c      = Jd + (n-1)*S;                                                 
  result = Contract(a,c);
end