function Equal_area_projection(a1_ang, a2_ang, a3_ang)
%  Equal-area projection
     
%     compute r for equal-area projection, both hemispheres will be plotted
%     a1
      [~,a1in]  = find(a1_ang(2,:)<=(0.5*pi));
      [~,a1out] = find(a1_ang(2,:)>(0.5*pi));
      r1(a1in)  = sqrt(2) * sin(a1_ang(2,a1in)./2);
      r1(a1out) = sqrt(2) * cos(a1_ang(2,a1out)./2);
%     a2
      [~,a2in]  = find(a2_ang(2,:)<=(0.5*pi));
      [~,a2out] = find(a2_ang(2,:)>(0.5*pi));
      r2(a2in)  = sqrt(2) * sin(a2_ang(2,a2in)./2);
      r2(a2out) = sqrt(2) * cos(a2_ang(2,a2out)./2);
%     a3       
      [~,a3in]  = find(a3_ang(2,:)<=(0.5*pi));
      [~,a3out] = find(a3_ang(2,:)>(0.5*pi));
      r3(a3in)  = sqrt(2) * sin(a3_ang(2,a3in)./2);
      r3(a3out) = sqrt(2) * cos(a3_ang(2,a3out)./2);
   
%     equal-area projections of a1, a2, a3   
%     a1      
      subplot(1,3,1);
      t = 0 : .01 : 2 * pi;
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
%     phi<=pi/2, plot red dots 
      polar(a1_ang(1,a1in),r1(a1in),'.r')
%     phi>pi/2, plot red dots 
      polar(a1_ang(1,a1out),r1(a1out),'.r')
      hold off
      title('a1')
      
%     a2   
      subplot(1,3,2);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a2_ang(1,a2in),r2(a2in),'.r')
      polar(a2_ang(1,a2out),r2(a2out),'.r')
      hold off
      title('a2')
      
%     a3   
      subplot(1,3,3);
      P = polar(t, ones(size(t)));
      set(P, 'Visible', 'off')
      hold on
      polar(a3_ang(1,a3in),r3(a3in),'.r')
      polar(a3_ang(1,a3out),r3(a3out),'.r')
      hold off
      title('a3')
end
