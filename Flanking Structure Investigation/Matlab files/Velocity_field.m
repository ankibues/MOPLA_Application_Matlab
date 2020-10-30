function [vel_f] = Velocity_field(ep,YY,ZZ,Vel_field)

% this is the function to calculate Velocity field at the current position
% of the marker points.

%Input: ep : 2*n matrix , coordinates of 'n' points 
%       YY and ZZ are the coordinates of the mesh
%       Vel_field , is the velocity field for the mesh; Obviously this will
%       be generated at every time step, and should have a separate
%       function.

[~,n]= size(ep);
vel_f = zeros(2,1,n);
 for i=1:n
    dist= ((ep(1,i)-ZZ).^2 + (ep(2,i)-YY).^2).^(0.5);
    minimum = min(min(dist));
    [xx,yy]=find(dist==minimum);
    vel_f(:,:,i)= Vel_field(:,:,xx,yy);
 end

end
