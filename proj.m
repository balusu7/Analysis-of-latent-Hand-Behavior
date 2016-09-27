clear;
F = dir('20*');
X=importdata(F(1).name,' ');
for ii = 2:length(F)
%import data from file
X=vertcat(X,importdata(F(ii).name,' '));
end
[U,S,normX,Z,W,Xmean] = ipca(X);
plot( Z(:,1),Z(:,2),'*');
hold on;
%testing with mat PCA
%[Uorg,Score,Eig]=pca(X);
%gmm on Z
AIC = zeros(1,4);
gm = cell(1,4);
%options = statset('MaxIter',10000);
for k = 1:3
    gm{k} = fitgmdist(Z,k);
    AIC(k)= gm{k}.AIC;
end
%finding min gmm required
[minAIC,numComponents] = min(AIC);
 gmm=fitgmdist(Z,4);%,'Replicates',10);%,'Start','plus');%,'RegularizationValue',0.1);%,'Options',options,'CovarianceType','diagonal');
ezcontour(@(x,y)pdf(gmm,[x y]),[-40 60],[-40 60])

% taking point from distribution
%mu=gmm.mu(1,1);
%sigma=gmm.Sigma(1,1);
%newZ =  pdf('Normal',mu+1,mu,sigma);
obj = gmdistribution(gmm.mu,gmm.Sigma,gmm.ComponentProportion);
Y = random(obj,1000);
scatter(Y(:,1),Y(:,2),10,'*r');

P = posterior(gmm,Y);
M = max(P,[],2);

%lets us trace back for one point
newZ=Y(1,:);
newX=newZ*W';
newX=newX+Xmean;
 newX=int8(newX);
oldZ=Z(1,:);
oldX=oldZ*W';
oldX=oldX+Xmean;

hold on;








