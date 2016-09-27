clear;
F = dir('20*');
X=importdata(F(1).name,' ');
for ii = 2:length(F)
%import data from file
X=vertcat(X,importdata(F(ii).name,' '));
end
[U,S,normX,Z,W,Xmean] = ipca(X);
%plot( Z(:,1),Z(:,2),'*');
X=Z;

% Set 'm' to the number of data points.
m = size(X, 1);

k = 3;  % The number of clusters.
n = 2;  % The vector lengths.

% Randomly select k data points to serve as the initial means.
indeces = randperm(m);
mu = X(indeces(1:k), :);
imu=mu;
sigma = [];

% Use the overal covariance of the dataset as the initial variance for each cluster.
for (j = 1 : k)
    sigma{j} = cov(X);
end

% Assign equal prior probabilities to each cluster.
p = ones(1, k) * (1 / k);
% Expectation Maximization

% Matrix to hold the probability that each data point belongs to each cluster.
% One row per data point, one column per cluster.
W = zeros(m, k);

% Loop until convergence.
for (iter = 1:100)
    % STEP 3a: Expectation
    %
    % Calculate the probability for each data point for each distribution.
    
    % Matrix to hold the pdf value for each every data point for every cluster.
    % One row per data point, one column per cluster.
    pdf = zeros(m, k);
    
    % For each cluster...
    for (j = 1 : k)
        
        % Evaluate the Gaussian for all data points for cluster 'j'.
        pdf(:, j) = mvnpdf(X, mu(j, :), sigma{j});%gaussianND(X, mu(j, :), sigma{j});
    end
    
    % Multiply each pdf value by the prior probability for cluster.
    %    pdf  [m  x  k]
    %    p  [1  x  k]   
    %  pdf_w  [m  x  k]
    Wpdf = bsxfun(@times, pdf, p);
    
    % Divide the weighted probabilities by the sum of weighted probabilities for each cluster.
    %   sum(pdf_w, 2) -- sum over the clusters.
    W = bsxfun(@rdivide, Wpdf, sum(Wpdf, 2));
    % STEP 3b: Maximization
    % Calculate the probability for each data point for each distribution.

    % Store the previous means.
    prevMu = mu;    
    
    % For each of the clusters...
    for (j = 1 : k)
    
        % Calculate the prior probability for cluster 'j'.
        p(j) = mean(W(:, j));
        
        % Calculate the new mean for cluster 'j' by taking the weighted
        % average of all data points.
        
        %mu(j, :) = weightedAverage(W(:, j), X);
        weights=W(:, j); values=X;
        %function [ val ] = weightedAverage(weights, values)
        val = weights' * values;

         % Divide by the sum of the weights.
          val = val ./ sum(weights, 1);
         
         mu(j, :)=val;
         
         
        % Calculate the covariance matrix for cluster 'j' by taking the 
        % weighted average of the covariance for each training example. 
        
        sigma_k = zeros(n, n);
        
        % Subtract the cluster mean from all data points.
        Xm = bsxfun(@minus, X, mu(j, :));
        
        % Calculate the contribution of each training example to the covariance matrix.
        for (i = 1 : m)
            sigma_k = sigma_k + (W(i, j) .* (Xm(i, :)' * Xm(i, :)));
        end
        
        % Divide by the sum of weights.
        sigma{j} = sigma_k ./ sum(W(:, j));
    end
    mu
    % Check for convergence.
    if (mu == prevMu)
        mu
        break
    end
            
% End of Expectation Maximization    
end
% STEP 4: Plot the data points and their estimated pdfs.
% 
% % Display a scatter plot of the two distributions.
% figure(2);
% hold off;
% plot(X1(:, 1), X1(:, 2), 'bo');
% hold on;
% plot(X2(:, 1), X2(:, 2), 'ro');
% 
% set(gcf,'color','white') % White background for the figure.
% 
% plot(mu1(1), mu1(2), 'kx');
% plot(mu2(1), mu2(2), 'kx');
% 
% % First, create a [10,000 x 2] matrix 'gridX' of coordinates representing
% % the input values over the grid.
% gridSize = 100;
% u = linspace(-6, 6, gridSize);
% [A B] = meshgrid(u, u);
% gridX = [A(:), B(:)];
% 
% % Calculate the Gaussian response for every value in the grid.
% z1 = gaussianND(gridX, mu(1, :), sigma{1});
% %z2 = gaussianND(gridX, mu(2, :), sigma{2});
% 
% % Reshape the responses back into a 2D grid to be plotted with contour.
% Z1 = reshape(z1, gridSize, gridSize);
% %Z2 = reshape(z2, gridSize, gridSize);
% 
% % Plot the contour lines to show the pdf over the data.
% [C, h] = contour(u, u, Z1);
% %[C, h] = contour(u, u, Z2);
% axis([-6 6 -6 6])
hold off;
plot( Z(:,1),Z(:,2),'.');
hold on;
% First, create a [10,000 x 2] matrix 'gridX' of coordinates representing
% the input values over the grid.
gridSize = 100;
u = linspace(-100, 100, gridSize);
[A,B] = meshgrid(u, u);
gridX = [A(:), B(:)];

% Calculate the Gaussian response for every value in the grid.
z1 = gaussianND(gridX, mu(1, :), sigma{1});
z2 = gaussianND(gridX, mu(2, :), sigma{2});
z3 = gaussianND(gridX, mu(3, :), sigma{3});

% Reshape the responses back into a 2D grid to be plotted with contour.
Z1 = reshape(z1, gridSize, gridSize);
Z2 = reshape(z2, gridSize, gridSize);
Z3 = reshape(z3, gridSize, gridSize);
% Plot the contour lines to show the pdf over the data.
[C, h] = contour(u, u, Z1,'k');
[C, h] = contour(u, u, Z2,'b');
[C, h] = contour(u, u, Z3,'r');
axis([-100 100 -100 100])

title('Original Data and Estimated PDFs');