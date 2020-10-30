function Ex= Extfield(a,EP,invS,u1,d1)
% This function calculates the exterior strain rate field at a given point
% from an inclusion, in an isotropic matrix

% a : length of axis
% EP : coordinate of the exterior point
% invS : inverse of S, S is  eshelby tensor
% u1 : interior perturbation of strain rate field i.e. dE1-d1
% d1 is macroscale strain rate field in inclusion coordinate system

% Exterior Field for 1st inclusion
        %--------------------------------------------------------------------------
        % G
        
        G1point = Ex_Gtensor(a,EP);
        
        % LMABDA_Ex 
        LAMBDA_Ex   = zeros(3,3);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex(ii,j) = -1/3.*(G1point(ii,j,1,1)+ G1point(ii,j,2,2)+ G1point(ii,j,3,3));
                LAMBDA_Ex(j,ii) = LAMBDA_Ex(ii,j);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G1point(1,1,iii,iii)+ G1point(2,2,iii,iii)+ G1point(3,3,iii,iii));
            LAMBDA_Ex(iii,iii) = squeeze(LAMBDA_Ex(iii,iii))+t;
        end
        
        % S_Ex & PI_Ex
        S_Ex   = zeros(3,3,3,3);
        PI_Ex  = zeros(3,3,3,3);
        delt   = eye(3);
        for ii=1:3
            for j=ii:3
                for k=1:3
                    for l=k:3
                        %S_Ex
                        S_Ex(ii,j,k,l)  = squeeze(G1point(ii,j,k,l))+ delt(k,l).*...
                                            squeeze(LAMBDA_Ex(ii,j));
                        S_Ex(j,ii,k,l)  = S_Ex(ii,j,k,l);
                        S_Ex(j,ii,l,k)  = S_Ex(ii,j,k,l);
                        S_Ex(ii,j,l,k)  = S_Ex(ii,j,k,l);
                        %PI_Ex
                        PI_Ex(ii,j,k,l) = 1/2.*(delt(j,k).*LAMBDA_Ex(ii,l) +...
                                            delt(j,l).*LAMBDA_Ex(ii,k) - delt(ii,k)...
                                            .*LAMBDA_Ex(j,l) - delt(ii,l).*...
                                            LAMBDA_Ex(j,k));
                        PI_Ex(ii,j,l,k) = PI_Ex(ii,j,k,l);
                        PI_Ex(j,ii,k,l) = -PI_Ex(ii,j,k,l);
                        PI_Ex(j,ii,l,k) = -PI_Ex(ii,j,k,l);
                    end
                end
            end
        end
        S_Ex(1,1,1,1) = -(S_Ex(1,1,2,2)+S_Ex(1,1,3,3));
        S_Ex(2,2,2,2) = -(S_Ex(2,2,1,1)+S_Ex(2,2,3,3));
        S_Ex(3,3,3,3) = -(S_Ex(3,3,1,1)+S_Ex(3,3,2,2));
  
        % Exterior field (e_Ex) at the centroid of 1st
        % inclusion due to 2nd inclusion
                
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex, invS);
            v2          = Multiply(v1, u1);
            Ex = v2 + d1;
end
