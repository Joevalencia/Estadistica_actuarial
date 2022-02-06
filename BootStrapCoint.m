function [r p alpha beta A Omega LL epsilons] = BootStrapCoint(kappas,alpha,beta,order,TestStat,Sig)

[~,m] = size(kappas);

T = m-1;
Dkap = kappas(:,2:m)-kappas(:,1:T);
if order == 0
    kappas = [kappas; ones(1,m)];
    X = [];
elseif order == 1
    kappas = [kappas; 1:m];
    X = ones(1,T);
elseif order == 2
    kappas = [kappas;(1:m).^2];
    X = [ones(1,T); 1:T];
else
    kappas = [kappas; (1:m).^3];
    X = [ones(1,T); 1:T; (1:T).^2];
end

[N,~] = size(kappas);

A = [];

StrapNo = 5000;
r = 0;
Continue = 1;

while Continue == 1
    
    epsilons = zeros(N-1,T);
    if order > 0
        if r == 0
            A = Dkap*X'*inv(X*X');
            for t=2:m
                epsilons(:,t-1) = kappas(1:(N-1),t) - (kappas(1:(N-1),t-1) + A*X(:,t-1));
            end
        else        
            A = (Dkap - alpha(:,1:r)*beta(:,1:r)'*kappas(:,1:T))*X'*inv(X*X');
            for t=2:m
                epsilons(:,t-1) = kappas(1:(N-1),t) - (kappas(1:(N-1),t-1) + A*X(:,t-1) + alpha(:,1:r)*beta(:,1:r)'*kappas(:,t-1));
            end
        end
    else
        if r == 0
            epsilons = Dkap;
        else
            for t=2:m
                epsilons(:,t-1) = kappas(1:(N-1),t) - (kappas(1:(N-1),t-1) + alpha(:,1:r)*beta(:,1:r)'*kappas(:,t-1));
            end
        end
    end
    Omega = cov(epsilons');
    C = chol(Omega);
    
    kappas_strap = zeros(N-1,m,StrapNo);
    TestStat_Strap = zeros(length(TestStat),StrapNo);
    for i = 1:StrapNo
        kappas_strap(:,:,i) = kappas(1:(N-1),:);
        for t=2:m
            if order > 0
                if r == 0
                    kappas_strap(:,t,i) = kappas_strap(:,t-1,i) + A*X(:,t-1) + C*randn(N-1,1);
                else
                    kappas_strap(:,t,i) = kappas_strap(:,t-1,i) + A*X(:,t-1) + alpha(:,1:r)*beta(:,1:r)'*[kappas_strap(:,t-1,i);kappas(N,t-1)] + C*randn(N-1,1);
                end
            else
                if r == 0
                    kappas_strap(:,t,i) = kappas_strap(:,t-1,i) + C*randn(N-1,1);
                else
                    kappas_strap(:,t,i) = kappas_strap(:,t-1,i) + alpha(:,1:r)*beta(:,1:r)'*[kappas_strap(:,t-1,i);kappas(N,t-1)] + C*randn(N-1,1);
                end
            end
        end
    [~,~,~,TestStat_Strap(:,i)] = Coint_3(kappas_strap(:,:,i),order);
    end
    p(r+1) = sum(TestStat_Strap(r+1,:) >= TestStat(r+1))/StrapNo;
    if p(r+1) > Sig
        if r == 0
            alpha = zeros(N-1,1);
            beta = zeros(N,1);
        else
            alpha = alpha(:,1:r);
            beta = beta(:,1:r);
        end
        LL = -(T/2)*log(det(Omega)) - trace((epsilons'/Omega)*epsilons);
        Continue = 0;
    else
        r = r + 1;
    end
end