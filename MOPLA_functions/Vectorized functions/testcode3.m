%testcode3 : checking two Av result
N=3;
I=1;
J=2;
[Jd, ~, ~, ~] = FourIdentity();
C= 2*10*Jd;
q   = zeros(3,3,N);
[A, ang] = RandAANG(10,N);
for i= 1:N
        
       q(:,:,i)= Q(ang(:,i));
      
end
   
Cm= Transform2(C,q,N);
[Alp,Bet,ww] = GaussGGQ(2);
theta = reshape(Alp,1,[]);
 phi   = reshape(Bet,1,[]);
 xi = [cos(theta).*sin(phi); sin(theta).*sin(phi); cos(phi)];

tic
Avinvvv1 = Avinv(Cm,xi);
toc

tic
Avinvvv2 = Avinv2(Cm,xi);
toc

function Avinvv = Avinv(Cm,xi)
  [~,~,~,~,n]=size(Cm);
  [~,N]=size(xi);
Aij= zeros(3,3,N);
Avinvv=zeros(4,4,n,N);
%Avinvvv=zeros(n,N);

for r=1:n
   for i=1:3
       for j=1:3
             %%%%%%% 
          Aij(i,j,:)= Cm(i,1,j,1,r)*xi(1,:).*xi(1,:) + Cm(i,1,j,2,r)*xi(1,:).*xi(2,:) +Cm(i,1,j,3,r)*xi(1,:).*xi(3,:)...
                    + Cm(i,2,j,1,r)*xi(2,:).*xi(1,:) + Cm(i,2,j,2,r)*xi(2,:).*xi(2,:) +Cm(i,2,j,3,r)*xi(2,:).*xi(3,:)...
                    + Cm(i,3,j,1,r)*xi(3,:).*xi(1,:) + Cm(i,3,j,2,r)*xi(3,:).*xi(2,:) +Cm(i,3,j,3,r)*xi(3,:).*xi(3,:);
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
  %{ 
 %%%%%-------inverse of Av matrix......this thing is taking time...try c code with this.  
  for rr= 1:N           
      Avinvv(:,:,rr) = inv(Av(:,:,rr));
  end
   %}
  Avinvv(:,:,r,:)= Av;
 
end

end



function Av = Avinv2(Cm,xi)
  [~,N]=size(xi);
  [~,~,~,~,n]= size(Cm);
Aij= zeros(3,3,n,N);
%Avinvv=zeros(4,4,N);
%Avinvvv2=zeros(n,N);

   for i=1:3
       for j=1:3
             %%%%%%% discuss this part
          Aij(i,j,:,:)= repmat(squeeze(Cm(i,1,j,1,:)),1,N).*repmat(xi(1,:).*xi(1,:),n,1) + repmat(squeeze(Cm(i,1,j,2,:)),1,N).*repmat(xi(1,:).*xi(2,:),n,1) +repmat(squeeze(Cm(i,1,j,3,:)),1,N).*repmat(xi(1,:).*xi(3,:),n,1)...
                    + repmat(squeeze(Cm(i,2,j,1,:)),1,N).*repmat(xi(2,:).*xi(1,:),n,1) + repmat(squeeze(Cm(i,2,j,2,:)),1,N).*repmat(xi(2,:).*xi(2,:),n,1) +repmat(squeeze(Cm(i,2,j,3,:)),1,N).*repmat(xi(2,:).*xi(3,:),n,1)...
                    + repmat(squeeze(Cm(i,3,j,1,:)),1,N).*repmat(xi(3,:).*xi(1,:),n,1) + repmat(squeeze(Cm(i,3,j,2,:)),1,N).*repmat(xi(3,:).*xi(2,:),n,1) +repmat(squeeze(Cm(i,3,j,3,:)),1,N).*repmat(xi(3,:).*xi(3,:),n,1);
       end
   end
      
 %%%%%--------- defining the Av matrix.
   Av(1,1,:,:)= Aij(1,1,:,:);
   Av(1,2,:,:)= Aij(1,2,:,:);
   Av(1,3,:,:)= Aij(1,3,:,:);
   Av(1,4,:,:)= repmat(xi(1,:),n,1);
   Av(2,1,:,:)= Aij(2,1,:,:);
   Av(2,2,:,:)= Aij(2,2,:,:);
   Av(2,3,:,:)= Aij(2,3,:,:);
   Av(2,4,:,:)= repmat(xi(2,:),n,1);
   Av(3,1,:,:)= Aij(3,1,:,:);
   Av(3,2,:,:)= Aij(3,2,:,:);
   Av(3,3,:,:)= Aij(3,3,:,:);
   Av(3,4,:,:)= repmat(xi(3,:),n,1);
   Av(4,1,:,:)= repmat(xi(1,:),n,1);
   Av(4,2,:,:)= repmat(xi(2,:),n,1);
   Av(4,3,:,:)= repmat(xi(3,:),n,1);
   Av(4,4,:,:)= 0;
   
 %{  
 for r=1:n  
 %%%%%-------inverse of Av matrix......this thing is taking time...try c code with this.  
  for rr= 1:N           
      Avinvv(:,:,rr) = inv(Av(:,:,rr));
  end
  
  Avinvvv2(r,:)= squeeze(Avinvv(I,J,:))';

 end
   %}
end
