clear;
hold off;
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
CovnX=cov(normX);
Y=X;%normX;
m=2;
[lambda,psi,T,stats,F] = factoran(X,m);
l=lambda*lambda';%covariance mattrix estimated
expcov=l + diag(psi);%expeted covariance
Xdata=X;
X=F;
Z=F;
siz = size(X, 1);
k = 3;  % The number of clusters.
n = 2;  % no of dimentions

% Randomly select k data points to serve as the initial means.
s = rng;
r = randi(siz,1,k)
%Assign first k random points to mu intialisation
mu = X(r, :);
imu=mu;
sigma = [];
%sigma is cell the inner contain nxn elements simliar to cov(x)
%variance across dimen and symmetric matrix

% Use the dataset as the initial variance for all clusters.
for (j = 1 : k)
    sigma{j} = cov(X);
end
isigma=sigma;
% Assign equal prior probabilities to each cluster.
p = ones(1, k);
p=p/k;
% Expectation Maximization


W = zeros(siz, k);%probability of each data point in k clusters 

% Loop until convergence.
for (iter = 1:100)
    %  Expectation
    
    pdf = zeros(siz, k);%for each iteration pdf is zero for all points in k clusters
    
    % For each cluster...
    for (j = 1 : k)
        Mu=mu(j, :);
        Sig=sigma{j};
        % Evaluate the Gaussian for all data points for cluster 'j'.
        pdf(:, j) = mvnpdf(X, Mu,Sig );%pdf of every x
    end
    
    % Multiply each pdf value by the prior probability for cluster.
    Wpdf = bsxfun(@times, pdf, p);
    
    % Divide wpdf by the sum of wpdf for all cluster.
   
    totalprobX=sum(Wpdf, 2);%row wise addition
    
    % final probabiltiy given x by suing bayes rule.
    W = bsxfun(@rdivide, Wpdf, totalprobX);
   
   
    % Calculate the probability for each data point for each distribution.

    % the previous means.
    prevMu = mu;    
    
    %  Maximization
    % For each of k clusters
    for (j = 1 : k)
    
        % the prior probability for cluster j.
        p(j) = mean(W(:, j));
        
        % the updated mean for cluster j by taking the average of corrsponding data points.
        
        weight=W(:, j);
        
        values=X;
        
        val = weight' * values;
        SumW=sum(weight);
         % Divide by the sum of the weights.
        val = val ./ SumW;%val is vector and element ny elemnt operation
         
         mu(j, :)=val;
         
         
        % Calculate the covariance matrix for cluster 'j' by taking the 
        % weighted average of the covariance for each training example. 
        
        SigmaK = zeros(n, n);
        
        % Subtract the cluster mean from all data points.
        Xm = bsxfun(@minus, X, mu(j, :));
        
        %  each training example to Sigma.
        for (i = 1 : siz)
            SigmaK = SigmaK + (W(i, j) .* (Xm(i, :)' * Xm(i, :)));
        end
        
        % Divide by the sum of weights to get fianl sigma for each cluster.
        sigma{j} = SigmaK ./ sum(W(:, j));
    end
    mu
    % Check for convergence.
    if (mu == prevMu)
        mu
        break
    end
               
end


%  Ploting  data  

plot( Z(:,1),Z(:,2),'.');
hold on;
% grid of size 100x100
gridSize = 100;
u = linspace(-2, 2, gridSize);
[A,B] = meshgrid(u, u);
gridX = [A(:), B(:)];

% Calculate the Gaussian response for every value in the grid.
z1 = mvnpdf(gridX, mu(1, :), sigma{1});
z2 = mvnpdf(gridX, mu(2, :), sigma{2});
z3 = mvnpdf(gridX, mu(3, :), sigma{3});

% Reshape from surface to 2-d plot contours ( concentric ellipse like )
Z1 = reshape(z1, gridSize, gridSize);
Z2 = reshape(z2, gridSize, gridSize);
Z3 = reshape(z3, gridSize, gridSize);
% Plot the contour lines to show the pdf over the data.
[C, h] = contour(u, u, Z1,'k');
[C, h] = contour(u, u, Z2,'b');
[C, h] = contour(u, u, Z3,'r');
axis([-2 2 -2 2])

% Sigma cell to 3 dim matrix
for d=1:k
    Sigma(:,:,d)=sigma{d};
end

obj = gmdistribution(mu,Sigma,p);
Y = random(obj,1000);
tscatter(Y(:,1),Y(:,2),10,'*r');
for ix=1:22
y(ix,1)=normrnd(0,psi(ix,1));
end
newX=F(1,:)*lambda'+psi';
newX=newX+Xmean;
diff=newX-Xdata(1,:);

%walk along axes
Yaxis= [-2:0.05:0];%0-80
Xaxis=zeros(size(Yaxis,2),1)';%(180+4*Yaxis)/-3;%0
Caxis = vertcat(Xaxis,Yaxis);
Caxis=Caxis';
Xsynaxi=Caxis*lambda';
Xsynaxis=zeros(size(Xaxis,2),22);
for i=1:size(Xaxis,2)
Xsynaxis(i,:)=Xsynaxi(i,:)+Xmean+psi';
end
%plot( Xaxis,Yaxis,'g*');
Xsynaxis=round(Xsynaxis);
dlmwrite('ppcawalkaxisY.txt',Xsynaxis,'delimiter','\t');
%walk along x axes
Xaxis= [-2:0.05:1.5];%0-80
Yaxis=zeros(size(Xaxis,2),1)';%(180+4*Yaxis)/-3;%0
Caxis = vertcat(Xaxis,Yaxis);
Caxis=Caxis';
Xsynaxi=Caxis*lambda';
Xsynaxis=zeros(size(Xaxis,2),22);
for i=1:size(Xaxis,2)
Xsynaxis(i,:)=Xsynaxi(i,:)+Xmean+psi';
end
plot( Xaxis,Yaxis,'g*');
Xsynaxis=round(Xsynaxis);
dlmwrite('ppcawalkaxisX.txt',Xsynaxis,'delimiter','\t');
