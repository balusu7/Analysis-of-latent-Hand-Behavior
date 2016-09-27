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
Y=X;%normX;
k=2;
[coeff,score,pcvar,mu,v,s] = ppca(Y,k);

y=normrnd(zeros(1,22),v,1,22);
newX=score(1,:)*coeff';
newX=newX+Xmean+y;
diff=newX-X(1,:);

plot(score(:,1),score(:,2),'.');
hold on;


Yaxis= [60:-2:-10];%0-80
Xaxis=zeros(size(Yaxis,2),1)';%(180+4*Yaxis)/-3;%0
Caxis = vertcat(Xaxis,Yaxis);
Caxis=Caxis';
%plot( Xaxis,Yaxis,'g*');
y=normrnd(zeros(1,22),v,1,22);
Xsynaxi=Caxis*coeff';
Xsynaxis=zeros(size(Xaxis,2),22);
for i=1:size(Xaxis,2)
Xsynaxis(i,:)=Xsynaxi(i,:)+Xmean+y;
end
plot( Xaxis,Yaxis,'g*');
Xsynaxis=round(Xsynaxis);
dlmwrite('newppcawalkaxisY.txt',Xsynaxis,'delimiter','\t');

Xaxis= [-80:2:50];%0-80
Yaxis=zeros(size(Xaxis,2),1)';%(180+4*Yaxis)/-3;%0
Caxis = vertcat(Xaxis,Yaxis);
Caxis=Caxis';
plot( Xaxis,Yaxis,'g*');
y=normrnd(zeros(1,22),v,1,22);
Xsynaxi=Caxis*coeff';
Xsynaxis=zeros(size(Xaxis,2),22);
for i=1:size(Xaxis,2)
Xsynaxis(i,:)=Xsynaxi(i,:)+Xmean+y;
end
%plot( Xaxis,Yaxis,'g*');
Xsynaxis=round(Xsynaxis);
dlmwrite('newppcawalkaxisX.txt',Xsynaxis,'delimiter','\t');
