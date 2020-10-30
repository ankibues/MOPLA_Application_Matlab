%testcode2  %%%% testing the Contract_vectorized and Contract_vectorized2 function %%%%
%%%% Latter can be used when all 'n' 4th order tensors to be multiplied by a
%%%% single 4th order tensor.
%%%% sample comparison with normal 'Contract' function. 

%%%% Elapsed time is 0.159954 seconds.
%%%% Elapsed time is 0.010096 seconds.
%%%% Elapsed time is 0.158095 seconds.
%%%% Elapsed time is 0.008371 seconds.


%%%% reducing the time by an order. would be useful while MOPLA operation!
 
 
%{
[Jd, ~, ~, ~] = FourIdentity();
Cm= 2*10*Jd;
a= sort(rand(3,1),'descend');
a= repmat(a,1,100);

T = SnpIn_all(a,Cm);

T1= T(:,:,:,:,1);
T2= T(:,:,:,:,2);
T3= T(:,:,:,:,3);


A= rand(4,4,100,20);
Avinvv= zeros(4,4,100,20);
tic
Ai= Ainv(A,100,20);
toc

tic
for r=1:100
 %%%%%-------inverse of Av matrix......this thing is taking time...try c code with this.  
   for rr= 1:20           
      Avinvv(:,:,r,rr) = inv(A(:,:,r,rr));
   end
 end

 toc
 %}

[Jd, ~, ~, ~] = FourIdentity();
cc= rand(3,3,3,3,200);
T= rand(3,3,3,3,200);
tic
z= zeros(3,3,3,3,200);
zzzz= zeros(3,3,3,3,200);
for r= 1:200
       z(:,:,:,:,r)= Contract(cc(:,:,:,:,r),T(:,:,:,:,r));
end
toc

tic
zz= Contract_vectorized(cc,T);
toc

tic 
for r= 1:200
       zzzz(:,:,:,:,r)= Contract(Jd,T(:,:,:,:,r));
end
toc

tic
zzz= Contract_vectorized2(Jd,T);
toc

