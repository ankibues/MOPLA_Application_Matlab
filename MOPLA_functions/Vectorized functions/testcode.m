%testcode  %%%% testing the Transform2 function %%%% vectorized equivalent
%of 'Transform'.
%%%% sample comparison with normal 'Transform' function. 
%%%% Elapsed time is 0.651465 seconds.
%%%% Elapsed time is 0.218307 seconds.


%%%% reducing the time by 1/3rd amount. would be useful while MOPLA operation!
%{
N= 200;
a1=10;

[a, ang] = RandAANG(a1,N);
q   = zeros(3,3,N);
for i=1:N
    q(:,:,i)= Q(ang(:,i));
         
end

X= rand(3,3,3,3);
 
ss = zeros(3,3,3,3,N);

tic
  for i=1:N
        
  ss(:,:,:,:,i)= Transform(X,q(:,:,i));
  end 
  toc
  
  
  
  tic
  s=Transform2(X,q,N);
  toc
  %}

%%%%% also testing contract1_vectorized
%%%%
% N= 100
%Elapsed time is 0.000959 seconds.
%Elapsed time is 0.002054 seconds.
%%Elapsed time is 0.000945 seconds.
%Elapsed time is 0.000143 seconds.
% as it can be seen, result is quite varying depending on the matrices to
% be multiplied
%{
N= 100;
a= rand(3,3);
b= rand(3,3,N);
k=zeros(1,N);
tic
for i = 1:N
    k(:,i) = contract1(b(:,:,i),a);
end
toc
tic
kk= contract1_vectorized(b,a);
toc
%}

%%%% test for multinv function ! it is faster than usual approach 
N=200;
M= rand(5,5,N);
tic
X = multinv(M);
toc
xx= zeros(5,5,N);
tic
for i=1:N
    xx(:,:,i) = eye(5)/M(:,:,i);
end
toc
