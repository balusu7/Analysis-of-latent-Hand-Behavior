clear;
F = dir('20*');
X=importdata(F(1).name,' ');
for ii = 2:length(F)
%import data from file
X=vertcat(X,importdata(F(ii).name,' '));
end
[U,S,normX,Z,Wpca,Xmean] = ipca(X);
%plot( Z(:,1),Z(:,2),'*');
S=S*S;% square of singular values to get eigen values
X=Z;
sum=zeros(1,22);
sum(1,1)=S(1,1);
for i = 2:22
    sum(1,i)=sum(1,i-1)+S(i,i);
end
S=round(S);
sum=round(sum);
var=sum./sum(1,22);
    