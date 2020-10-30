function theta = freeslip(a,w,d,S)
% Wd.m
%
% the Ellipsoid vorticity referred to the frame tracking its semi-axes
%--------------------------------------------------------------------------
 
  theta = zeros(3,3);
  
  for i = 1:3
      for j = 1:3
          if i==j
              theta(i,j)=0;
          else
              k = d(i,j)*(((a(i)^2-a(j)^2)*2*S(i,j,i,j)/(a(i)^2+a(j)^2))-((a(i)^2+a(j)^2)/(a(i)^2-a(j)^2)))/(1-2*S(i,j,i,j));
              theta(i,j) = k + w(i,j);
          end
      end
  end
  
end
