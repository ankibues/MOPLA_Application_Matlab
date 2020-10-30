function Vel_field= VelGradtoVelfield(a,YY,ZZ,nn,ll,LL)

Vel_field= zeros(2,1,nn,nn);
       
        
        % -------getting the velocity field using integration approach.
       % first quadrant of mesh
        for i= (nn+1)/2 :(nn-1)      % this corresponds to movement along Y-direction of mesh
            
            for j= (nn+1)/2 :(nn-1)  % this corresponds to movement along X-direction of mesh
            ind= (YY(i,j+1)/a(2))^2 + (ZZ(i,j+1)/a(3))^2;
            if ind<1
                Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j+1)-ZZ(i,j)); 
                Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j+1)-ZZ(i,j));
            end
            if ind>1
                Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + LL(2,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
                Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + LL(3,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
            end
            end
            ind= (YY(i+1,(nn+1)/2)/a(2))^2 + (ZZ(i+1,(nn+1)/2)/a(3))^2;
            if ind<1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            if ind>1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            
        end
       
        % second quadrant of mesh
            for i= (nn+1)/2 :-1:2      % this corresponds to movement along Y-direction of mesh
                  for j= (nn+1)/2 :(nn-1)  % this corresponds to movement along X-direction of mesh
                     ind= (YY(i,j+1)/a(2))^2 + (ZZ(i,j+1)/a(3))^2;
                     if ind<1
                        Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j+1)-ZZ(i,j)); 
                        Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j+1)-ZZ(i,j));
                     end
                     if ind>1
                        Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + LL(2,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
                        Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + LL(3,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
                     end
                  end
                  ind= (YY(i-1,(nn+1)/2)/a(2))^2 + (ZZ(i-1,(nn+1)/2)/a(3))^2;
                  if ind<1
                     Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                     Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                  end
                  if ind>1
                     Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                     Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                  end
            
            end
        
% third quadrant of mesh
        for i= (nn+1)/2 :-1:2      % this corresponds to movement along Y-direction of mesh
            
            for j= (nn+1)/2 :-1:2  % this corresponds to movement along X-direction of mesh
            ind= (YY(i,j-1)/a(2))^2 + (ZZ(i,j-1)/a(3))^2;
            if ind<1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j-1)-ZZ(i,j)); 
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j-1)-ZZ(i,j));
            end
            if ind>1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + LL(2,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + LL(3,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
            end
            end
            ind= (YY(i-1,(nn+1)/2)/a(2))^2 + (ZZ(i-1,(nn+1)/2)/a(3))^2;
            if ind<1
                Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            if ind>1
                Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            
       end
     %fourth quadrant of mesh
     
        for i= (nn+1)/2 :(nn-1)      % this corresponds to movement along Y-direction of mesh
            
            for j= (nn+1)/2 :-1:2  % this corresponds to movement along X-direction of mesh
            ind= (YY(i,j-1)/a(2))^2 + (ZZ(i,j-1)/a(3))^2;
            if ind<1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j-1)-ZZ(i,j)); 
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j-1)-ZZ(i,j));
            end
            if ind>1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + LL(2,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + LL(3,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
            end
            end
            ind= (YY(i+1,(nn+1)/2)/a(2))^2 + (ZZ(i+1,(nn+1)/2)/a(3))^2;
            if ind<1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            if ind>1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            
        end
end

        