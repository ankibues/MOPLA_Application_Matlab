function ed = Ed(Nm,Nc,r,s,d,epsilonII,b,Jd)

% Iterative calculation in the ellipsoid, starting with current strain rate
% invariant "REF", until the new strain rate invariant converges. Both the
% strain rate tensor "de" and the viscosity ratio (final yy) and the strain
% rate invariant are output.
%**********This procedure is to consider power law
%materials************incompressible isotropic power law material**********

% Input: Nm,Nc are power-law stress exponents for the matrix and ellipsoid, respectively;
%        r is the viscosity ratio(at matrix strain rate), 
%        epsilonII is the strain rate invariant at which ellipsoid
%        viscosity is defined;
%        s is the Eshelby tensor S, a 4th order tensor;
%        d is the strain rate tensor(for the matrix as I understand) defined in the ellipsoidal system at
%        current state, 3*3 matrix;(updated 1/2/2017 AB)
%        
% 28/03/2015


  yy = r;
  REF = epsilonII;
  for jj = 1:10000
      de       = Multiply(fdE(Nm,yy,s,b,Jd),d); % eq. 16 of Jiang 2012b
      epsilonI = inva(de);
      zz       = (epsilonI/REF)^((1-Nc)/Nc)*yy;   %section 5.4 Jiang2012b
      alpha    = abs(zz/yy - 1);
      yy       = zz; % new r value
      REF      = epsilonI; % new strain rate invariant value
      if alpha < 0.001
          break
      end
  end
  nn  = [yy;REF;0];
  ed = [de,nn];           % new strain rate tensor, new viscosity ratio, new strain rate invariant
end

function result = fdE(n,rr,S,b,Jd)
  
  a      = Inverse(Jd+((n*rr)-1)*S,b);    % eq. 16 of Jiang 2012b(deformable ellipsoid in power law viscous materials)
  c      = Jd + (n-1)*S;                    % but with Jd instead of Js( that is done for incompressible materials)                              
  result = Contract(a,c);
end

function rr = inva(x)
 rr = (0.5*contract1(x,x))^0.5;
end