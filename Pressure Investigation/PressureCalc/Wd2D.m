function wp = Wd2D(a,m,d)
% Ellipsoid vorticity referred to the frame tracking its semi-axes

% Input: a is two semi-axes, 2*1 matrix
%        m is the W in the ellipsoid coordinate, 2*2 matrix;
%        d is the new strain rate tensor, 2*2 matrix;

% Output: wp updated vorticity, 2*2 matrix;

 
  wp = zeros(2,2);          % wp is simply the shear spin part of rotation of ellipsoid.
  
  for i = 1:2
      for j = 1:2
          if i == j         % find the reason why?
              wp(i,j) = 0;
          elseif a(i) == a(j)
              wp(i,j) = m(i,j);
          else
              wp(i,j) = p(a(i),a(j))*d(i,j);
          end
      end
  end
  
end

function r = p(x,y)
  
r = (x^2+y^2)/(x^2-y^2);

end
