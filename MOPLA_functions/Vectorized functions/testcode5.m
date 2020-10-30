% test for multiply_vectorized function
%for N=200 
%Elapsed time is 0.017499 seconds.
%Elapsed time is 0.001459 seconds.
%
%for N=100 
%Elapsed time is 0.009377 seconds.
%Elapsed time is 0.002089 seconds.

% Not bad ! Quite faster than the original one !

N=200;
a= rand(3,3);
b= rand(3,3,3,3,N);
m= zeros(3,3,N);
tic
for i= 1:N
    m(:,:,i) = Multiply(b(:,:,:,:,i),a);
end
toc

tic
M= Multiply_vectorized(b,a);
toc

% test for multiply_vectorized1 function
%N=100
%Elapsed time is 0.009579 seconds.
%Elapsed time is 0.004738 seconds.

%N=200
%Elapsed time is 0.017851 seconds.
%Elapsed time is 0.001596 seconds.