%create ER adjacency matrix
clear all

N=1000;   % number of nodes = size of the matrix
c=4;      % mean degree
p=c/N;    % link probability


%matrix generation
J=rand(N)<p; %generate logical matrix of 1(condition ok) and 0(not ok);
J=triu(J,1); %select triangle above the 1st diagonal-> main diag=0
J=J+J';      %add J to its transpose J'

%uncomment the following lines to get a weighted matrix with Gaussin weights
% K=triu(randn(N),1);
% K=K+K';
% J=J.*K;   %element-wise product

%J is an ER adjacency matrix.

%Uncomment the following to 
%check that the degree distribution of J is Poisson.

% %Expected deg distribution: Poisson with parameter c.
% pk=poisspdf(0:15,c);   %adjust the degree range if needed
% 
% %Actual degree distribution of J
% deg=sum(J~=0,2); 
% 
% %Compare pk with histogram of deg.
% histogram(deg,'Normalization','probability');
% hold on
% plot(0:15,pk,'*')
