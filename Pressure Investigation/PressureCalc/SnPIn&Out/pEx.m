function p = pEx(sq,a,x,Alp,Bet,ww)
 % Green interaction tensor and Eshelby Tensors inside and outside the
 % ellipsoid in an isotropic medium, according to Jiang's mathcad sheet
 % SnPI(a,x)
 
 % here, only Green Interaction tensor for Pressure considered for
 % different quadratures. for other, see SnpIn.m file.
 
 % Input: a is three semi-axes of the inclusion(a1>a2>a3), 3*1 matrix;
 %        x is the coordinates of the point that we considered, 3*1 matrix;
 %        sq is the sphere quard method used
 %           choices:apgq, ggq, ggqs, glesh; if the global methods(ggq,ggqs,glesh) 
 %           are chosen, there is a number of points(gp) in varargin
 
 % Output: 
 %         p is a 2nd order tensor that is related to pressure.
 
 % 29/03/2015
 
   xx = [a,x]; 
   p = zeros(3,3);
   T = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k= 1:3
                for l= 1:3
                    if strcmp(sq,'apgq')
                        p(i,j) = APGQ('pgex',[i,j],xx,0,2*pi,0,pi,10^(-4));
                    elseif strcmp(sq,'glesh')
                        p(i,j) = GLeSH('pgin',[i,j],a,Alp,Bet,ww);   % just considering Internal pressure
                    elseif strcmp(sq,'ggq')
                        T(i,j,k,l) = GGQ('fgin',[i,j,k,l],a,0,2*pi,0,pi,Alp,Bet,ww); % just considering internal stress
                    else
                        fprintf('need to choose a method (APGQ, GLeSH, GGQ)');
                    end
                end    
            end
        end
    end
  
end

