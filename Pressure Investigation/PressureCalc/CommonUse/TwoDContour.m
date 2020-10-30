function TwoDContour(P)
 % contour plotting for a 2D contour
 
 
 % Input: P is a n*3 matrix, column1 is interested value at n ponits, column2&3 are 
 %                           correspounding x,y coordinates of each point.
 %                           x,y are positive integers; y is the row number
 %                           of the point, and x is the column number.
 %        
 
 % Output: a plot
 
   
   p    = P(:,1);
   xx   = P(:,2);
   yy   = P(:,3);

   xmax = max(xx);
   ymax = max(yy);
 
   P   = accumarray([yy,xx],p,[ymax,xmax]);
 
   figure
 
   pcolor(1:xmax,1:ymax,P);
   shading interp;
   axis image
   axis off
   title('Pressure')
   colorbar('vert')
  
end