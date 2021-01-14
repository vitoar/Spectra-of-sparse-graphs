%Spectrum of a single instance
clear all


N=1000; % number of nodes=matrix size
c=4;    % mean degree
p=c/N;  % link probability



epsi=10^-3;          % Epsilon (regulariser)
delta_lambda=0.01;   %lambda scan
%It can be chosen smaller than epsi if a super high resolution is needed.
lambdamax=3; %Need to be adjusted according to the matrix. 
%These values are ok for rescaled ER matrices.

lambdavalues=linspace(0,lambdamax,lambdamax/delta_lambda);
lambda_size=length(lambdavalues);
sample_size=1;  
%for a single instance choose n=1. To get some(useful)averaging effect,choose n>1.
rho=zeros(lambda_size,sample_size);

for n=1:sample_size 
    
    
    %Matrix generation: here for instance we consider {0,1}-ER matrices: 
    %Choose a different generator according to the needs. 
    %You can also generate matrices and load them in this script 
    %before entering this loop. 
    %Python's networkx package contains useful adjcency matrices generators.
    
    J=rand(N)<p; %generate logical matrix of 1(condition ok) and 0(not ok);
    J=triu(J,1); %select triangle above the 1st diagonal-> main diag=0
    J=J+J';
    
    
    %rescaling matrix elements
    J=J./sqrt(c);
   
    
    %degree vector
    deg=sum(J~=0,2); 
     
 
    %neighbours matrix
    NB=zeros(N,max(deg)+2);
    for i=1:N
        NB(i,1)=i;
        NB(i,2)=deg(i);
        if NB(i,2)~=0 
            v=find(J(i,:));
            for j=1:NB(i,2)
                NB(i,2+j)=v(j);
            end
        end
    
    end
    
    
    
    
 for s=1:lambda_size
    lambda=lambdavalues(s);


    %initilize cavity inverse variances with Re[w]>0
    w=zeros(N);

    for i=1:N
        if NB(i,2)~=0
            for j=1:NB(i,2)
            w(i,NB(i,2+j))=10+rand+1i*randn;
            end
        end
    end




    %Self consistency solution for cavity marginals
    w_vect_old=w(:);
    w_vect_old(w_vect_old==0)=[];
    t_w=0;
    epsilon_w=10^(-10);
    w_ratio=2;

    %cavity inverse variance routine
    while w_ratio>epsilon_w

        for i=1:N
        if NB(i,2)~=0
            for j=1:NB(i,2)
                sum_w=0;
                for l=1:NB(NB(i,2+j),2)
                    if NB(NB(i,2+j),2+l)~=i
                        sum_w=sum_w+J(NB(i,2+j),NB(NB(i,2+j),2+l))*J(NB(i,2+j),NB(NB(i,2+j),2+l))/w(NB(i,2+j),NB(NB(i,2+j),2+l));
                    
                    end
                end

            w(i,NB(i,2+j))=1i*lambda+epsi+sum_w;
      
            end
        end
        end
  
% Compute ratio R to assess convergence
    w_vect_curr=w(:);
    w_vect_curr(w_vect_curr==0)=[];
    w_ratio=norm((w_vect_curr-w_vect_old),1)/norm(w_vect_old,1);
    w_vect_old=w_vect_curr;

    t_w=t_w+1;
    end

%test for non negative real part
test=all(all(real(w)>=0));
if test==0
    disp('step')
    disp(t_w)
    return
end


% marginal inverse variances
W=zeros(N,1);


 for i=1:N
    sumW=0;
    if NB(i,2)~=0  
        for j=1:NB(i,2)
            sumW=sumW+J(i,NB(i,2+j))^2/w(i,NB(i,2+j));
        end
    end
    W(i)=1i*lambda+epsi+sumW;
 end
 

 
    rho(s,n)=mean(real(1./W))/pi;
 
 
 end
 
 
end

ave_rho=mean(rho,2);  %This is half of the asd.

%Mirror the asd (even function of lambda) to negative values of lambda.
neg_rho=flip(ave_rho);
neg_rho(end)=[];
neg_lambdavalues=linspace(-lambdamax,0,lambdamax/delta_lambda);
neg_lambdavalues(end)=[];
lv=[neg_lambdavalues';lambdavalues'];
tot_rho=[neg_rho;ave_rho];   %Full asd

plot(lv,tot_rho) %Plot full asd
