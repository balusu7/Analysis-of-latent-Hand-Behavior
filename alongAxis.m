clear;
F = dir('20*');
X=importdata(F(1).name,' ');
for ii = 2:length(F)
%import data from file
X=vertcat(X,importdata(F(ii).name,' '));
end
[U,S,normX,Z,Wpca,Xmean] = ipca(X);
%plot( Z(:,1),Z(:,2),'*');
X=Z;

% Set 'm' to the number of data points.
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
    %mu
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
u = linspace(-100, 100, gridSize);
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
axis([-100 100 -100 100])

% Sigma cell to 3 dim matrix
for d=1:k
    Sigma(:,:,d)=sigma{d};
end

obj = gmdistribution(mu,Sigma,p);
Y = random(obj,1000);
%scatter(Y(:,1),Y(:,2),10,'*r');

P = W;
M = max(P,[],2);

%lets us trace back for one point
newZ=Y(1,:);
newX=newZ*Wpca';
newX=newX+Xmean;

Xsyn=Y*Wpca';
for i=1:1000
Xsyn(i,:)=Xsyn(i,:)+Xmean;
end
% walking along training data Xsyndata
Xsyndat=X*Wpca';
Xsyndata=zeros(100,22);
for i=1:100
Xsyndata(i,:)=Xsyndat(i,:)+Xmean;
end
Xsyndata=round(Xsyndata);
dlmwrite('pcawalkdata.txt',Xsyndata,'delimiter','\t');
%walking along axes Xsynaxis
Xaxis=[0:-1:-60];%0-80
Yaxis=(180+4*Xaxis)/-3;%0zeros(size(Xaxis,2),1)';%
Caxis = vertcat(Xaxis,Yaxis);
Caxis=Caxis';
Xsynaxi=Caxis*Wpca';
Xsynaxis=zeros(size(Xaxis,2),22);
for i=1:size(Xaxis,2)
Xsynaxis(i,:)=Xsynaxi(i,:)+Xmean;
end
Xsynaxis=round(Xsynaxis);
dlmwrite('pcawalkaxisrandom.txt',Xsynaxis,'delimiter','\t');
plot( Xaxis,Yaxis,'g*');
%posterior of new data 
pdf1 = zeros(size(Xaxis,2), k);
 for (j = 1 : k)
        Mu=mu(j, :);
        Sig=sigma{j};
        pdf1(:, j) = mvnpdf(Caxis, Mu,Sig );
 end
    Wpdf = bsxfun(@times, pdf1, p);
    totalprobX=sum(Wpdf, 2);
    W = bsxfun(@rdivide, Wpdf, totalprobX);
   PnewX=W;
  PnewX = max(PnewX,[],2);

  if min(PnewX)<0.8
      fprintf('probability decreased along this axis');
  elseif fprintf('probability is good along this axis');
  end