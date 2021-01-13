%create ER adjacency matrix
clear all

N=1000;   % number of nodes = size of the matrix
c=4;      % mean degree
p=c/N;    % link probability


%matrix generation
J=rand(N)<p; %generate logical matrix of 1(condition ok) and 0(not ok);
J=triu(J,1); %select triangle above the 1st diagonal-> main diag=0
J=J+J';      %add J to its transpose J'

%uncomment to get a weighted matrix with Gaussin weights
K=triu(randn(N)./sqrt(c),1);
K=K+K';
J=J.*K;   %element-wise product

