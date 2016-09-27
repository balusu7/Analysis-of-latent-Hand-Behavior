F = dir('20*');
X=importdata(F(1).name,' ');
Xmean=mean(X);
for i=1:size(X,1)
    
    normX(i,:)=X(i,:)-Xmean;
    
end
CovnX=cov(normX);
[U,S,V]=svd(CovnX);
d=1;%d is name of row of index of row
Z=normX*U;
newX=Z(d,:)*U';
newX=newX+Xmean;
diff=newX-X(d,:);