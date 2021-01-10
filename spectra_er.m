%population dynamics per spettri

clear all

lambdamax=3.5;
lambdavalues=linspace(0,lambdamax,lambdamax*100);
lambda_size=length(lambdavalues);
rho=zeros(lambda_size,1);
error=zeros(lambda_size,1);
epsi=10^-6;


c=4;
N=4000;

kmin=0;
kmax=13;
%probability values for any given k, p(k)
p=zeros(kmax-kmin+1,1);
for k=kmin:kmax
    p(k+1)=poisspdf(k,c)/poisscdf(kmax,c); 
end


deg_seq=kmin:kmax;
k_avg=sum(deg_seq'.*p);

%r(k)=kp(k)/<k>
r=zeros(kmax,1);
for k=1:kmax
    r(k)=k*p(k+1)/k_avg;
end

for s=1:lambda_size
    s
    lambda=lambdavalues(s);
    lambdae=lambda-1i*10^-300;
    
    %Equilibration*******************************************************
    E=500;
    re_m_w=zeros(E,1);
    im_m_w=zeros(E,1);
    w=10+rand(N,1)+1i*randn(N,1); %initial condition with Re(w)>0
    %mi serve una nuova popolazione per ogni nuovo valore di lambda
    for e=1:E
        
        %degree sequence generation
        deg1=[];
        for i=1:kmax
            deg1=[deg1; repmat(i,floor(N*r(i)),1)];
        end
        while length(deg1)<N
            deg1=[deg1; deg1(randi(length(deg1)))];
        end
        deg1=deg1(randperm(N));
        
        %equilibration sweep
        for i=1:N
            k=deg1(i);
            ind=randperm(N,k-1);
            J=randn/sqrt(c);
            %J=1;
            w(i)=J^2/(1i*lambdae+sum(w(ind)));
        end
        re_m_w(e)=mean(real(w));
        im_m_w(e)=mean(imag(w));
    end
    
    %Measurement*********************************************************
    %A sweep before actual measure
    
    M=10^4; %number of meas sweeps
    n=10^7/M; %number of samples per meas sweep
    x=zeros(M*n,1);
    a=zeros(M*n,1);
    b=zeros(M*n,1);
    
    %generate degree sequence according to p(k)
        deg2=[];
        for i=kmin:kmax
            deg2=[deg2; repmat(i,floor(M*n*p(i+1)),1)];
        end
        while length(deg2)<M*n
            deg2=[deg2; deg2(randi(length(deg2)))];
        end
        deg2=deg2(randperm(M*n));
    
    for m=1:M
        
    
    
        %degree sequence generation
        deg1=[];
        for i=1:kmax
            deg1=[deg1; repmat(i,floor(N*r(i)),1)];
        end
        while length(deg1)<N
            deg1=[deg1; deg1(randi(length(deg1)))];
        end
        deg1=deg1(randperm(N));
        
        %sweep S(lambda)
        for i=1:N
            k=deg1(i);
            ind=randperm(N,k-1);
            J=randn/sqrt(c);
            %J=1;
            w(i)=J^2/(1i*lambdae+sum(w(ind)));
        end
    
        %measurement

        for j=1:n
            k=deg2((m-1)*n+j);
            ind=randperm(N,k);
            x((m-1)*n+j)=sum(w(ind));
            a((m-1)*n+j)=real(x((m-1)*n+j));
            b((m-1)*n+j)=imag(x((m-1)*n+j));
        end

        
    end
    
    rho(s)=mean((a+epsi)./(pi*((a+epsi).^2.+(b+lambda).^2)));
    error(s)=std((a+epsi)./(pi*((a+epsi).^2.+(b+lambda).^2)));
end

neg_rho=flip(rho);
neg_rho(end)=[];
neg_lambdavalues=linspace(-lambdamax,0,lambdamax*100);
neg_lambdavalues(end)=[];
lv=[neg_lambdavalues';lambdavalues'];
tot_rho=[neg_rho;rho];
plot(lv,tot_rho)

save('ER_c4_Gaussian4000_epsi106_mln_k13.mat')