%Population Dynamics algorithm for the average spectral density (asd) 
%of ER matrices

clear all

epsi=10^-20;   %epsilon to be used for the asd in the sampling procedure
lambdamax=3.5;  %adjust it according to needs
delta_lambda=0.01;  %Ok to get the full asd. 
%If you want to resolve a small interval in lambda with high res, 
%then delta_lambda must be a fraction of epsi.
%E.g. delta_lambda=0.1*espi;
lambdavalues=linspace(0,lambdamax,lambdamax/delta_lambda);
lambda_size=length(lambdavalues);
rho=zeros(lambda_size,1);


c=4;      %mean degree
N=1000;   %population size

kmin=0;    %min degree
kmax=18;   %for practical purpose, choose a max degree

%degree distribution p(k)
p=zeros(kmax-kmin+1,1);
for k=kmin:kmax
    p(k+1)=poisspdf(k,c)/poisscdf(kmax,c); 
end


deg_seq=kmin:kmax;
k_avg=sum(deg_seq'.*p);  %check for the actual mean degree

%r(k)=kp(k)/<k> 
r=zeros(kmax,1);
for k=1:kmax
    r(k)=k*p(k+1)/k_avg;
end


for s=1:lambda_size
    
    lambda=lambdavalues(s);
    lambdae=lambda-1i*10^-300;
    
    %Equilibration sweeps*****************
    E=500;
    re_m_w=zeros(E,1);
    im_m_w=zeros(E,1);
    w=10+rand(N,1)+1i*randn(N,1); %initial condition with Re(w)>0
    %a new population needed for any new value of lambda
    
    %degree sequence generation (once for all equilibration sweeps)
    deg1=[];
    for i=1:kmax
        deg1=[deg1; repmat(i,floor(N*E*r(i)),1)];
    end
    while length(deg1)<N*E
        deg1=[deg1; deg1(randi(length(deg1)))];
    end
    deg1=deg1(randperm(N*E));
    
    
    for e=1:E
        % main part of a single equilibration sweep
        for i=1:N
            k=deg1((e-1)*N+i);
            ind=randperm(N,k-1);
            J=randn/sqrt(c);  %gaussian weights
            %J=1/sqrt(c);     %unweighted case
            w(i)=J^2/(1i*lambdae+sum(w(ind)));
        end
        re_m_w(e)=mean(real(w)); %to check convergence. E=500 seems enough here.
        im_m_w(e)=mean(imag(w)); %as above
        
    end
    
    %Measurement Sweeps *******************
    
    
    M=10^4; %number of meas sweeps
    n=10^7/M; %number of samples per meas sweep
    x=zeros(M*n,1);
    a=zeros(M*n,1);
    b=zeros(M*n,1);
    
    %generate degree sequence according to p(k)(once for all!)
        deg2=[];
        for i=kmin:kmax
            deg2=[deg2; repmat(i,floor(M*n*p(i+1)),1)];
        end
        while length(deg2)<M*n
            deg2=[deg2; deg2(randi(length(deg2)))];
        end
        deg2=deg2(randperm(M*n));
        
        
        %degree sequence generation for meas sweeps (once for all!)
        deg1=[];
        for i=1:kmax
            deg1=[deg1; repmat(i,floor(M*N*r(i)),1)];
        end
        while length(deg1)<N*M
            deg1=[deg1; deg1(randi(length(deg1)))];
        end
        deg1=deg1(randperm(M*N));
    
    for m=1:M
        
    
    
       
        %A sweep  S(lambda) before actual measure
        for i=1:N
            k=deg1((m-1)*N+i);
            ind=randperm(N,k-1);
            J=randn/sqrt(c);
            %J=1/sqrt(c);
            w(i)=J^2/(1i*lambdae+sum(w(ind)));
        end
    
        %measurement sweep

        for j=1:n
            k=deg2((m-1)*n+j);
            ind=randperm(N,k);
            x((m-1)*n+j)=sum(w(ind));
            a((m-1)*n+j)=real(x((m-1)*n+j));
            b((m-1)*n+j)=imag(x((m-1)*n+j));
        end

        
    end
    
    
    rho(s)=mean((a+epsi)./(pi*((a+epsi).^2.+(b+lambda).^2)));
    
end


%mirror the asd (even function of lambda) to negative values of lambda.
neg_rho=flip(rho);
neg_rho(end)=[];
neg_lambdavalues=linspace(-lambdamax,0,lambdamax*100);
neg_lambdavalues(end)=[];
lv=[neg_lambdavalues';lambdavalues'];
tot_rho=[neg_rho;rho];

plot(lv,tot_rho) %plot asd
