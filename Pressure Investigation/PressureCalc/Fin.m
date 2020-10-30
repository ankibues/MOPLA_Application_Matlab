function T = Fin(a,Alp,Bet,ww)
 % Green interaction tensor for  and Eshelby Tensors inside and outside the
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
 
   
   
   T = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k= 1:3
                for l= 1:3
                    T(i,j,k,l) = GGQ('fgin',[i,j,k,l],a,0,2*pi,0,pi,Alp,Bet,ww); % just considering internal stress
                end    
            end
        end
    end
  
end

