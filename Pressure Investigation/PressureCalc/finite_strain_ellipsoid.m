function [eval,evec]=finite_strain_ellipsoid(F)
%Input:Position gradient tensor
%Output:eval: three principal strains  eval=[e1;e2;e3], e1>e2>e3
%       evec: three principal strain axes  evec=[v1 v2 v3] 
%             vi and ei are paired.

% Polar decomposition  
  B     = F * F';
% Eigenvectors and corresponding eigenvalues of B
  [V,D] = eig(B);
% Calculate eigenvalues
  eval  = diag(D).^0.5;
% Sort the eigenvalues in descending order  
  [~,permutation] = sort(eval,'descend');
  eval            = eval(permutation);
% Rearrange the eigenvectors
  V            = V(:,permutation);
  V(:,3)       = cross(V(:,1),V(:,2));
  evec         = V;
end