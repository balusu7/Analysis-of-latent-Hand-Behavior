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

%[F,adj_var,cum_var] = sparsePCA(X, card, num_comp, num_runs, verbosity)
[F,adj_var,cum_var] = sparsePCA(normX, 15, 2, 1000, 1);
Y=normX*F;
newX=Y*F';
%newx=newX+Xmean;
