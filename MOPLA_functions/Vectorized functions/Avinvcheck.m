function A = Avinvcheck(I,J,Cm,xi)
  [~,N]=size(xi);
Aij= zeros(3,3,N);
Avinv=zeros(4,4,N);

        for i=1:3
            for j=1:3
                for m=1:3 %%%%%%% discuss this part....not sure how m and n decided.
                    for n=1:3
                Aij(i,j,:)= Cm(i,m,j,n).*xi(m,:).*xi(n,:);
                    end
                end
            end
        end
      
 %%%%%--------- defining the Av matrix.
   Av(1,1,:)= Aij(1,1,:);
   Av(1,2,:)= Aij(1,2,:);
   Av(1,3,:)= Aij(1,3,:);
   Av(1,4,:)= xi(1,:);
   Av(2,1,:)= Aij(2,1,:);
   Av(2,2,:)= Aij(2,2,:);
   Av(2,3,:)= Aij(2,3,:);
   Av(2,4,:)= xi(2,:);
   Av(3,1,:)= Aij(3,1,:);
   Av(3,2,:)= Aij(3,2,:);
   Av(3,3,:)= Aij(3,3,:);
   Av(3,4,:)= xi(3,:);
   Av(4,1,:)= xi(1,:);
   Av(4,2,:)= xi(2,:);
   Av(4,3,:)= xi(3,:);
   Av(4,4,:)= 0;
   
 %%%%%-------inverse of Av matrix  
  for r= 1:N
      Avinv(:,:,r) = inv(Av(:,:,r));
  end
  
  A= squeeze(Avinv(I,J,:))';
  
end