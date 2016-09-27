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
Z=normX';
[Zica A T mu] = myICA(Z,2);

newX=Zica(:,1)'*A;
newX=newX+mu'+Xmean;

