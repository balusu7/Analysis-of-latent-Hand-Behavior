clear;
F = dir('20*');
X=importdata(F(1).name,' ');
for ii = 2:length(F)
%import data from file
X=vertcat(X,importdata(F(ii).name,' '));
end
Xmean=mean(X);
%normalizing the data
for i=1:size(X,1)
    for j=1:size(X,2)
    normX(i,j)=X(i,j)-Xmean(j);
    end
end
% normX is zero mean
CovnX=cov(X);
Y=X;%normX;
m=2;
[lambda,psi,T,stats,F] = factoran(X,m);
y=normrnd(zeros(1,22),psi',1,22);
newX=F(2,:)*lambda';
newX=newX+Xmean+y;
diff=newX-X(2,:);

