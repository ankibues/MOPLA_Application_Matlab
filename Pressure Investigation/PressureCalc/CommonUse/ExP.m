function outp = ExP(x,y,h,u2,sq,Alp,Bet,ww)

% This function is to get pressure deviation at y 
% According to Jiang's mathcad sheet StateI(x,y) in EshalbyTensorIn&Out3.

% Input: x is the initial state of the system, a 3*6 matrix-- x(1:3,1:3) is
% the orientation of the inclusion(q0); x(1:3,4) is the three semi-axes of
% the inclusion; x(1,5) is the viscosity ratio (at matrix strain rate),
% x(2,5) is the strain rate invariant at which ellipsoid viscosity is
% defined, x(3,5) is 0; x(1:3,6) is the pressure deviation, which are all 0
% at initial state.
%        y is a point we consisdered, a 3*1 matrix.
%        h = Contract(2*Jd,invS);
%        u2 is the strain rate difference between clast and matrix;
%        sq is the sphere quard method used
%           choices:apgq, ggq, ggqs, glesh; if the global methods(ggq,ggqs,glesh) 
%           are chosen, there is a number of points(gp) in varargin

% Output

% Update:28/03/2015

  a     = x(:,4);
  
  pe    = pEx(sq,a,y,Alp,Bet,ww);
  p_out = R_Multiply(pe,h);
  
  outP  = contract1(p_out,u2)*[1;0;0];                                    
  outp  = outP(1);

end