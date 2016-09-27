clear;
F = dir('20*');
X=importdata(F(1).name,' ');
X=X';Z = kernelpca_tutorial(X,3);
%  Ploting  data  

plot( Z(:,1),Z(:,2),'b*');figure;plot( Z(:,2),Z(:,3),'b*');figure;plot( Z(:,3),Z(:,1),'b*');

hold on;
