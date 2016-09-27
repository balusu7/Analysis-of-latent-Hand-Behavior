function [U,S,normX,Z,W,Xmean] = ipca(X)

% caluculating mean for each dimention
Xmean=mean(X);
%normalizing the data
for i=1:size(X,1)
    for j=1:size(X,2)
    normX(i,j)=X(i,j)-Xmean(j);
    end
end
%Computing covariance matrix
CovnX=cov(normX);
%extracting eigen components using SVD U eigen vector dec order and 
%S eigen diagonal vector in dec order
[U,S,V]=svd(CovnX);
%checking orthogonolity of U (just testing) 
orthoU=U*U';
%testing
%U=U*-1;
%extracting 2 (low) dim from high 22
W=U(:,1:1:2); % w is 22x2
Z=normX*W; % 111x22 22x2 = 111x2









