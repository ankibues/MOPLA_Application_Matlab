function ellipsoidplot(a1,a2,a3)
%a1=10;
%a2=10;
%a3=1;
[x, y, z] = ellipsoid(0,0,0,a1,a2,a3,1000);
%figure
s= surf(x, y, z);
%set(s, 'FaceColor',[.83 .82 .78])


%ellipsoidplot(3,1.5,1)s.FaceColor = 'gray';

shading 'flat';
material 'metal';
lighting 'gouraud';
%axis tight;
camlight('local')
%rotate(s,[1 0 0],0)

end


