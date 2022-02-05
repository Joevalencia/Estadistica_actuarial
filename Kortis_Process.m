DataFit;

SimNo = 5000;%100;%10000;%2000;%
ProjPeriod = 25;
set_seed(1981);

WeightsEWg = ZeroCohorts(YearsEW,AgesEW,WeightsEWm,10);
WeightsUSg = ZeroCohorts(YearsUS,AgesUS,WeightsUSm,10);
muEWm = CreateMuUni(alphaEWm,kappasEWm,fxEWm,zeros(1,mEW),zeros(nEW,1),gammasEWm);
muUSm = CreateMuUni(alphaUSm,kappasUSm,fxUSm,zeros(1,mUS),zeros(nUS,1),gammasUSm);

pEW_2011 = exp(-muEWm(:,mEW));
PEW_2011 = ones(nEW+1,1);
pUS_2011 = exp(-muUSm(:,mEW));
PUS_2011 = ones(nEW+1,1);
for x=2:(nEW+1)
    PEW_2011(x) = PEW_2011(x-1)*pEW_2011(x-1);
    PUS_2011(x) = PUS_2011(x-1)*pUS_2011(x-1);
end
dEW_2011 = -diff(PEW_2011);
DEW_2011 = zeros(nEW,1);
dUS_2011 = -diff(PUS_2011);
DUS_2011 = zeros(nEW,1);
for x=1:nEW
    DEW_2011(x) = sum(dEW_2011(1:x));
    DUS_2011(x) = sum(dUS_2011(1:x));
end

DEW = [ones(mEW-1,1);flipud(DEW_2011)]';
DUS = [ones(mEW-1,1);flipud(DUS_2011)]';

%%

gammas_proj = zeros(2,cEW+ProjPeriod,SimNo);
gammas_projEW = zeros(SimNo,cEW+ProjPeriod);
gammas_projUS = zeros(SimNo,cEW+ProjPeriod);

DummyYears = [1919 1920 1921 1946 1947];
X_dummies = double((ones(5,1)*CohortsEW)==(DummyYears'*ones(1,cEW)))-double((ones(5,1)*CohortsEW)==((DummyYears-1)'*ones(1,cEW)));
mu_dummiesEW = gammasEWm(10:(cEW-9))*X_dummies(:,10:(cEW-9))'*inv(X_dummies(:,10:(cEW-9))*X_dummies(:,10:(cEW-9))');
gammasEWm = gammasEWm - mu_dummiesEW*X_dummies;
mu_dummiesUS = gammasUSm(10:(cEW-9))*X_dummies(:,10:(cEW-9))'*inv(X_dummies(:,10:(cEW-9))*X_dummies(:,10:(cEW-9))');
gammasUSm = gammasUSm - mu_dummiesUS*X_dummies;

f = @(Param) -GamLike(Param,gammasEWm(10:(cEW-9)),DEW(10:(cEW-9)),3);
ParamEW =[0.9 0.01 zeros(1,4)];
PEW = fminsearch(f, ParamEW);
rhoEW_gam = PEW(1);
sigmaEW_eps = PEW(2);
nuEW_0 = PEW(3:end);
[AlphaEWm KappasEWm GammasEWm] = ChangeGamID(alphaEWm,kappasEWm,gammasEWm,fxEWm,nuEW_0);

g = @(Param) -GamLike(Param,gammasUSm(10:(cEW-9)),DUS(10:(cEW-9)),3);
ParamUS =[0.9 0.01 zeros(1,4)];
PUS = fminsearch(g, ParamUS);
rhoUS_gam = PUS(1);
sigmaUS_eps = PUS(2);
nuUS_0 = PUS(3:end);
[AlphaUSm KappasUSm GammasUSm] = ChangeGamID(alphaUSm,kappasUSm,gammasUSm,fxUSm,nuUS_0);

%%

% Pairs = zeros(7,2);
% Pairs(:,1) = (1:7)';
% Pairs(:,2) = [1 2 5 6 3 7 4]';

% ProjPeriod = 25;

order = [3 2 1];
Sig = 0.05;

A = zeros(2,3,3);
Pi = zeros(2,3,3);
beta0 = zeros(3,3,3);
alpha0 = zeros(2,3,3);
beta = zeros(3,3,3);
alpha = zeros(2,3,3);
lambda = zeros(3,3);
TestStat = zeros(3,3);
r = zeros(1,3);
p = zeros(4,3);
Omega = zeros(2,2,3);
LL = zeros(1,3);
epsilons = zeros(2,mEW-1,3);

for i=1:3
    [beta0(:,:,i) alpha0(:,:,i) lambda(:,i) TestStat(:,i)] = Coint_3([KappasEWm(i,:);KappasUSm(i,:)],order(i));
    [r(i) p(1:(r(i)+1),i) alpha(:,1:max(r(i),1),i) beta(:,1:max(r(i),1),i) A(:,1:order(i),i) Omega(:,:,i) LL(i) epsilons(:,:,i)] = BootStrapCoint([KappasEWm(i,:);KappasUSm(i,:)],alpha0(:,:,i),beta0(:,:,i),order(i),TestStat(:,i),Sig);
    Pi(:,:,i) = alpha(:,:,i)*beta(:,:,i)';
end

%%
A_comb = zeros(10,3);
Pi_comb = zeros(10,14);
kappasEnd = [KappasEWm(:,mEW);KappasUSm(:,mEW);1;mEW;mEW^2;mEW.^3];
epsilons_comb = zeros(10,mEW-1);
for i=1:3
    A_comb(i,1:order(i)) = A(1,1:order(i),i);
    A_comb(i+4,1:order(i)) = A(2,1:order(i),i);
    Pi_comb(i,i) = Pi(1,1,i);
    Pi_comb(i,i+4) = Pi(1,2,i);
    Pi_comb(i+4,i) = Pi(2,1,i);
    Pi_comb(i+4,i+4) = Pi(2,2,i);
    Pi_comb(i,11+order(i)) = Pi(1,3,i);
    Pi_comb(i+4,11+order(i)) = Pi(2,3,i);
    epsilons_comb(i,:) = epsilons(1,:,i);
    epsilons_comb(i+4,:) = epsilons(2,:,i);
end

[A_comb(4,1),~,~,epsilons_comb(4,:)] = FitBrokenRW(YearsEW,KappasEWm(4,:),[]);
for i=4:6
    [A_comb(i+4,1),~,~,epsilons_comb(i+4,:)] = FitBrokenRW(YearsUS,KappasUSm(i,:),[]);
end

V = cov(epsilons_comb',0);
C = chol(V);
Z_stan = epsilons_comb;
for t=1:(mEW-1);
    Z_stan(:,t) = C\epsilons_comb(:,t);
end

%%

% mEW2009 = DeathsEWm(26:36,mEW-2)./ExposuresEWm(26:36,mEW-2);
% mUS2009 = DeathsUSm(6:16,mUS-2)./ExposuresUSm(6:16,mEW-2);
muEWm_proj = zeros(nEW,mEW+ProjPeriod,SimNo);
mEWm_proj = zeros(nEW,ProjPeriod,SimNo);
muUSm_proj = zeros(nUS,mUS+ProjPeriod,SimNo);
mUSm_proj = zeros(nUS,ProjPeriod,SimNo);

IndEW = zeros(SimNo,(mEW-8+ProjPeriod));
IndUS = zeros(SimNo,(mEW-8+ProjPeriod));
LDIVHist = zeros(SimNo,(mEW-8+ProjPeriod));

IndexEW = zeros(1,SimNo);
IndexUS = zeros(1,SimNo);
LDIV = zeros(1,SimNo);

kappas_proj = zeros(10,ProjPeriod,SimNo);
X_kap = [ones(1,ProjPeriod);(mEW+(1:ProjPeriod));(mEW+(1:ProjPeriod)).^2];

YearsX = YearsEW(end) + (1:ProjPeriod);
Year2009 = find(YearsX==2009);
Year2016 = find(YearsX==2016);
Years = YearsEW(8) + (1:(mEW-8+ProjPeriod));

mEWm = DeathsEWm./ExposuresEWm;
mUSm = DeathsUSm./ExposuresUSm;

for i=1:SimNo
    
    gammas_projEW(i,1:10) = GammasEWm(1:10);
    gammas_projUS(i,1:10) = GammasUSm(1:10);
    
    for y=10:(cEW-10)
        gammas_projEW(i,y) = DEW(y)*GammasEWm(y) + (1-DEW(y))*rhoEW_gam*gammas_projEW(i,y-1) + mu_dummiesEW*X_dummies(:,y) + sqrt((1-DEW(y)))*sigmaEW_eps*randn(1);
        gammas_projUS(i,y) = DUS(y)*GammasUSm(y) + (1-DUS(y))*rhoUS_gam*gammas_projUS(i,y-1) + mu_dummiesUS*X_dummies(:,y) + sqrt((1-DUS(y)))*sigmaUS_eps*randn(1);
    end

    for y=(cEW-9):(cEW+ProjPeriod)
        gammas_projEW(i,y) = rhoEW_gam*gammas_projEW(i,y-1) + sigmaEW_eps*randn(1);
        gammas_projUS(i,y) = rhoUS_gam*gammas_projUS(i,y-1) + sigmaUS_eps*randn(1);
    end
    
    for t=1:ProjPeriod
        Z = zeros(10,1);
        i_rand = floor(1+10*rand(1,10));
        t_rand = floor(1+(mEW-1)*rand(1,10));
        for j=1:10
            Z(j) = Z_stan(i_rand(j),t_rand(j));
        end
        if t==1
            kappas_proj(:,t,i) = kappasEnd(1:10,1) + A_comb*X_kap(:,t) + Pi_comb*kappasEnd + C*Z;
        else
            kappas_proj(:,t,i) = kappas_proj(:,t-1,i) + A_comb*X_kap(:,t) + Pi_comb*[kappas_proj(:,t-1,i);1;(mEW+t-1);(mEW+t-1)^2;(mEW+t-1)^3] + C*Z;
        end
    end
    
    muEWm_proj(:,:,i) = CreateMuUni(AlphaEWm,[KappasEWm kappas_proj(1:4,:,i)],fxEWm,zeros(1,mEW+ProjPeriod),zeros(nEW,1),gammas_projEW(i,:));
    muUSm_proj(:,:,i) = CreateMuUni(AlphaUSm,[KappasUSm kappas_proj(5:10,:,i)],fxUSm,zeros(1,mEW+ProjPeriod),zeros(nEW,1),gammas_projUS(i,:));
%     mEWm_proj(:,:,i) = ProjectDeaths(ExposuresEWm(:,mEW),DeathsEWm(:,mEW),muEWm_proj(:,(mEW+1):(mEW+ProjPeriod),i),ProjPeriod);
%     mUSm_proj(:,:,i) = ProjectDeaths(ExposuresUSm(:,mEW),DeathsUSm(:,mEW),muUSm_proj(:,(mEW+1):(mEW+ProjPeriod),i),ProjPeriod);
    mEWm_proj(:,:,i) = muEWm_proj(:,(mEW+1):(mEW+ProjPeriod),i);
    mUSm_proj(:,:,i) = muUSm_proj(:,(mEW+1):(mEW+ProjPeriod),i);
    
    [LDIVHist(i,:) IndEW(i,:) IndUS(i,:)] = CalcLDIV(AgesEW,[mEWm mEWm_proj(:,:,i)],[mUSm mUSm_proj(:,:,i)]);
    IndexEW(i) = IndEW(i,Years==2016);
    IndexUS(i) = IndUS(i,Years==2016);
    LDIV(i) = LDIVHist(i,Years==2016); 
    i
end

% LEDiff = LEEW - LEUS;

figure
hist(LDIV)

Levels = 0.034:0.001:0.039;
Exceed = zeros(1,6);
for i=1:6
    Exceed(i) = sum(LDIV >= Levels(i));
end

PRP = max(zeros(1,SimNo),min((LDIV-Levels(1))/(Levels(6)-Levels(1)),ones(1,SimNo)));
ExLoss = sum(PRP)/sum(PRP > zeros(1,SimNo));

%%
% IndEW = zeros(SimNo,(mEW-8+ProjPeriod));
% IndUS = zeros(SimNo,(mEW-8+ProjPeriod));
% LDIVHist = zeros(SimNo,(mEW-8+ProjPeriod));
% 
% mEWm = DeathsEWm./ExposuresEWm;
% mUSm = DeathsUSm./ExposuresUSm;
% ImpEW = zeros(11,(mEW-8));
% ImpUS = zeros(11,(mUS-8));
% for t=9:mEW
%     ImpEW(:,t-8) = 1 - (mEWm(26:36,t)./mEWm(26:36,t-8)).^(1/8);
%     ImpUS(:,t-8) = 1 - (mUSm(6:16,t)./mUSm(6:16,t-8)).^(1/8);
% end
% IndEW(:,1:(mEW-8)) = ones(SimNo,1)*mean(ImpEW,1);
% IndUS(:,1:(mEW-8)) = ones(SimNo,1)*mean(ImpUS,1);
% LDIVHist(:,1:(mEW-8)) = IndEW(:,1:(mEW-8)) - IndUS(:,1:(mEW-8));
% 
% ImpEWs = zeros(11,ProjPeriod,SimNo);
% ImpUSs = zeros(11,ProjPeriod,SimNo);
% for t=1:8
%     for i=1:SimNo
%         ImpEWs(:,t,i) = 1 - (mEWm_proj(26:36,t,i)./mEWm(26:36,(mEW-9+t))).^(1/8);
%         ImpUSs(:,t,i) = 1 - (mUSm_proj(6:16,t,i)./mUSm(6:16,(mUS-9+t))).^(1/8);
%     end
% end
% for t=9:ProjPeriod
%     ImpEWs(:,t,:) = 1 - (mEWm_proj(26:36,t,:)./mEWm_proj(26:36,t-8,:)).^(1/8);
%     ImpUSs(:,t,:) = 1 - (mUSm_proj(6:16,t,:)./mUSm_proj(6:16,t-8,:)).^(1/8);
% end
% for i=1:SimNo
%     IndEW(i,(mEW-7):(mEW-8+ProjPeriod)) = mean(ImpEWs(:,:,i),1);
%     IndUS(i,(mEW-7):(mEW-8+ProjPeriod)) = mean(ImpUSs(:,:,i),1);
%     LDIVHist(i,(mEW-7):(mEW-8+ProjPeriod)) = IndEW(i,(mEW-7):(mEW-8+ProjPeriod)) - IndUS(i,(mEW-7):(mEW-8+ProjPeriod));
% end
% Years = YearsEW(8) + (1:(mEW-8+ProjPeriod));
figure
FanChart([], [], Years, LDIVHist', 0.01:0.01:0.99);
hold
plot(Years,[0.039;0.034]*ones(1,length(Years)),'k--')
plot([2016 2016],[0.05 -0.03],'k:')

%%
ColoursEW = hsv(4);

figure
for i=1:4
    FanChart(YearsEW,KappasEWm(i,:)',YearsX,reshape(kappas_proj(i,:,:),ProjPeriod,SimNo),0.02:0.01:0.98,ColoursEW(i,:));
    hold
end

ColoursUS = hsv(6);
figure
for i=1:6
    FanChart(YearsUS,KappasUSm(i,:)',YearsX,reshape(kappas_proj(i+4,:,:),ProjPeriod,SimNo),0.02:0.01:0.98,ColoursUS(i,:));
    hold
end

CohortsX = [CohortsEW CohortsEW(end)+(1:ProjPeriod)];
figure
FanChart([],[],CohortsX,gammas_projEW', 0.02:0.01:0.98);
hold
plot(CohortsEW(1:(cEW-9)),GammasEWm(1:(cEW-9)) + mu_dummiesEW*X_dummies(:,1:(cEW-9)),'k:','linewidth',2)
hold

figure
FanChart([],[],CohortsX,gammas_projUS', 0.02:0.01:0.98);
hold
plot(CohortsUS(1:(cEW-9)),GammasUSm(1:(cEW-9)) + mu_dummiesUS*X_dummies(:,1:(cEW-9)),'k:','linewidth',2)
hold

% figure
% FanChart([],[],YearsX,LEEW',0.01:0.01:0.99,ColoursEW(1,:));
% hold
% FanChart([],[],YearsX,LEUS',0.01:0.01:0.99,ColoursEW(3,:));
% hold
% v = axis;
% v(3:4) = [25 40];
% axis(v);
% 
% figure
% FanChart([],[],YearsX,LEDiff',0.01:0.01:0.99);