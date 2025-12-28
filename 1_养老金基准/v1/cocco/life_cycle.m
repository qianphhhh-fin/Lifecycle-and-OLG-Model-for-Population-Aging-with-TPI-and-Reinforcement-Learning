%% Life-Cycle Consumption/Portfolio Choice Problem

clear all;
close all;
%% Variable Definitions

filename = string(zeros(80,1));

tb     = 20;
tr     = 66;
td     = 100;
tn     = td-tb+1;

na     = 51;
ncash  = 51;
n      = 5;
nc     = 21;

maxcash     = 200.0;
mincash     = 0.25;
aa          = -2.170042+2.700381;
b1          = 0.16818;
b2          = -0.0323371/10;
b3          = 0.0019704/100;
ret_fac     = 0.68212;
smay        = 0.1;
smav        = 0.1;
corr_v      = 0.0;
corr_y      = 0.0;
rho         = 10.0;
delta       = 0.97;
psi         = 0.5;
r           = 1.015;
mu          = 0.04;
sigr        = 0.2;

survprob    = zeros(tn-1,1);
delta2      = zeros(tn-1,1);
grid        = zeros(n,1);
weig        = zeros(n,1);
gret        = zeros(n,1);
ones_n_1    = ones(n,1);
grid2       = zeros(n,1);
yp          = zeros(n,n);
yh          = zeros(n,n);
nweig1      = zeros(n,n,n);
f_y         = zeros(tr-tb+1,1);
gy          = zeros(tr-tb,1);
gyp         = zeros(n,n,tn-1);
gcash       = zeros(ncash,1);
lgcash      = zeros(ncash,1);
ga          = zeros(na,1);
riskret     = zeros(na,n);
gc          = zeros(nc,1);
auxV        = zeros(na,nc);
vec_V       = zeros(na*nc,1);
secd        = zeros(ncash,1);
C           = zeros(ncash,tn);
V           = zeros(ncash,tn);
A           = ones(ncash,tn);

nsim        = 10000;
ones_nsim_1 = ones(nsim,1);
meanY       = zeros(tn,1);
meanC       = zeros(tn,1);
meanW       = zeros(tn,1);
meanA       = zeros(tn,1);
meanS       = zeros(tn,1);
meanB       = zeros(tn,1);
meanWY      = zeros(tn,1);
meanalpha   = zeros(tn,1);
meanGPY     = zeros(tn,1);
cGPY        = zeros(tn,1);
meanYs      = zeros(tn,1);
meanCs      = zeros(tn,1);
meanWs      = zeros(tn,1);
simPY       = zeros(tn,nsim);
simGPY      = zeros(tn,nsim);
simY        = zeros(tn,nsim);
simC        = zeros(tn,nsim);
simW        = zeros(tn,nsim);
simA        = zeros(tn,nsim);
simS        = zeros(tn,nsim);
simB        = zeros(tn,nsim);
simW_Y      = zeros(tn,nsim);
simR        = zeros(tn,nsim);
%% Approximation to Normal Distribution

grid(1,1) = -2.85697001387280;
grid(2,1) = -1.35562617997427;
grid(3,1) = 0.00000000000000;
grid(4,1) = 1.35562617997426;
grid(5,1) = 2.85697001387280;

weig(1,1) = 0.01125741132772;
weig(2,1) = 0.22207592200561;
weig(3,1) = 0.53333333333333;
weig(4,1) = 0.22207592200561;
weig(5,1) = 0.01125741132772;

%% Conditional Survival Probabilities

survprob(1,1)  = 0.99845;
survprob(2,1)  = 0.99839;
survprob(3,1)  = 0.99833;
survprob(4,1)  = 0.9983;
survprob(5,1)  = 0.99827;
survprob(6,1)  = 0.99826;
survprob(7,1)  = 0.99824;
survprob(8,1)  = 0.9982;
survprob(9,1)  = 0.99813;
survprob(10,1) = 0.99804;
survprob(11,1) = 0.99795;
survprob(12,1) = 0.99785;
survprob(13,1) = 0.99776;
survprob(14,1) = 0.99766;
survprob(15,1) = 0.99755;
survprob(16,1) = 0.99743;
survprob(17,1) = 0.9973;
survprob(18,1) = 0.99718;
survprob(19,1) = 0.99707;
survprob(20,1) = 0.99696;
survprob(21,1) = 0.99685;
survprob(22,1) = 0.99672;
survprob(23,1) = 0.99656;
survprob(24,1) = 0.99635;
survprob(25,1) = 0.9961;
survprob(26,1) = 0.99579;
survprob(27,1) = 0.99543;
survprob(28,1) = 0.99504;
survprob(29,1) = 0.99463;
survprob(30,1) = 0.9942;
survprob(31,1) = 0.9937;
survprob(32,1) = 0.99311;
survprob(33,1) = 0.99245;
survprob(34,1) = 0.99172;
survprob(35,1) = 0.99091;
survprob(36,1) = 0.99005;
survprob(37,1) = 0.98911;
survprob(38,1) = 0.98803;
survprob(39,1) = 0.9868;
survprob(40,1) = 0.98545;
survprob(41,1) = 0.98409;
survprob(42,1) = 0.9827;
survprob(43,1) = 0.98123;
survprob(44,1) = 0.97961;
survprob(45,1) = 0.97786;
survprob(46,1) = 0.97603;
survprob(47,1) = 0.97414;
survprob(48,1) = 0.97207;
survprob(49,1) = 0.9697;
survprob(50,1) = 0.96699;
survprob(51,1) = 0.96393;
survprob(52,1) = 0.96055;
survprob(53,1) = 0.9569;
survprob(54,1) = 0.9531;
survprob(55,1) = 0.94921;
survprob(56,1) = 0.94508;
survprob(57,1) = 0.94057;
survprob(58,1) = 0.9357;
survprob(59,1) = 0.93031;
survprob(60,1) = 0.92424;
survprob(61,1) = 0.91717;
survprob(62,1) = 0.90922;
survprob(63,1) = 0.90089;
survprob(64,1) = 0.89282;
survprob(65,1) = 0.88503;
survprob(66,1) = 0.87622;
survprob(67,1) = 0.86576;
survprob(68,1) = 0.8544;
survprob(69,1) = 0.8423;
survprob(70,1) = 0.82942;
survprob(71,1) = 0.8154;
survprob(72,1) = 0.80002;
survprob(73,1) = 0.78404;
survprob(74,1) = 0.76842;
survprob(75,1) = 0.75382;
survprob(76,1) = 0.73996;
survprob(77,1) = 0.72464;
survprob(78,1) = 0.71057;
survprob(79,1) = 0.6961;
survprob(80,1) = 0.6809;
%% Output Files - Names

for t=1:tn-1
    if t<10
        filename(t,1) = string('year0')+string(t);
    else
        filename(t,1) = string('year')+string(t);
    end
end
%% Additional Computations

for i1=1:n
    gret(i1,1) = r+mu+grid(i1,1)*sigr;
end

for i6=1:n
    for i7=1:n
        for i8=1:n
            nweig1(i6,i7,i8) = weig(i6,1)*weig(i7,1)*weig(i8,1);
        end
    end
end

theta = (1.0-rho)/(1.0-1.0/psi);
psi_1 = 1.0-1.0/psi;
psi_2 = 1.0/psi_1;
%% Grids for the State Variables and for Portfolio Rule

for i1=1:na
    ga(i1,1)=(na-i1)/(na-1.0);
end

for i5=1:na
    for i8=1:n
        riskret(i5,i8)=r*(1-ga(i5,1))+gret(i8,1)*ga(i5,1);
    end
end

l_maxcash = log(maxcash);
l_mincash = log(mincash);
stepcash = (l_maxcash-l_mincash)/(ncash-1);

for i1=1:ncash
    lgcash(i1,1)=l_mincash+(i1-1.0)*stepcash;
end
for i1=1:ncash
    gcash(i1,1)=exp(lgcash(i1,1));
end
%% Labor Income

for i1=1:n
    grid2(:,1) = grid(i1,1)*corr_y+grid(:,1).*ones_n_1(:,1)*(1-corr_y^2)^(0.5);
    yh(1:n,i1) = exp(grid2(:,1)*smay);
end

for i1=1:n
    grid2(:,1) = grid(i1,1)*corr_v+grid(:,1).*ones_n_1(:,1)*(1-corr_v^2)^(0.5);
    yp(:,i1) = grid2(:,1)*smav;
end

for i1=tb:tr
    f_y(i1-tb+1,1) = exp(aa+b1*i1+b2*i1^2+b3*i1^3);
end

for i1=tb:tr-1
    gy(i1-tb+1,1) = f_y(i1-tb+2,1)/f_y(i1-tb+1,1)-1.0;
    for i2=1:n
        gyp(:,i2,i1-tb+1) = exp(gy(i1-tb+1,1)*ones_n_1(:,1)+yp(:,i2));
    end
end

for i1=tr-tb+1:tn-1
    for i2=1:n
        gyp(:,i2,i1) = exp(0.0*ones_n_1(:,1));
    end
end
%% Terminal Period

for i1=1:ncash
    C(i1,tn) = gcash(i1,1);
end
A(:,tn) = 0.0;
for i1=1:ncash
    V(i1,tn) = C(i1,tn)*((1.0-delta)^(psi/(psi-1.0)));
end

secd(:,1) = f_spline(gcash(:,1),V(:,tn),ncash,1.0);
%% Retirement Periods

for i1 = 1:td-tr
    t = tn-i1;
    display(t)
    for i3=1:ncash
        if i3==1
            maxc = C(i3,t+1);
            minc = maxc/2;
        else
            minc = C(i3-1,t);
            if i3<10
                maxc = minc + (gcash(i3,1) - gcash(i3-1,1));
            else
                mpc = max((C(i3-1,t)-C(i3-9,t))/(gcash(i3-1,1) - gcash(i3-9,1)),0.1);
                maxc = minc + mpc*(gcash(i3,1) - gcash(i3-1,1));
            end
        end
        stepc=(maxc-minc)/(nc-1);
        for i9=1:nc
            gc(i9,1)=minc+(i9-1.0)*stepc;
        end
        for i4=1:nc
            u=(1.0-delta)*(gc(i4,1)^psi_1) ;
            sav = gcash(i3,1)-gc(i4,1);
            for i5=1:na
                auxVV=0.0;
                %           for i6=1:n
                %              for i7=1:n
                for i8=1:n
                    cash_1 = riskret(i5,i8)*sav+ret_fac;
                    cash_1 = max(min(cash_1,gcash(ncash,1)),gcash(1,1));
                    int_V  = f_sc_splint(gcash(:,1),V(:,t+1),secd(:,1),ncash,cash_1);
                    auxVV=auxVV+weig(i8,1)*survprob(t,1)*(int_V^(1.0-rho));
                end
                %              end
                %           end
                auxV(i5,i4) = (u+delta*(auxVV^(1.0/theta)))^psi_2;    %1/(1-1/psi)
            end
        end
        vec_V = reshape(auxV,[na*nc,1]);
        [V(i3,t),pt] = max(vec_V(:,1));
        aux2 = floor((pt(1)-1)/na);
        C(i3,t) = gc(aux2+1,1);
        A(i3,t) = ga(pt(1)-na*aux2,1);
    end

    secd(:,1) = f_spline(gcash(:,1),V(:,t),ncash,1.0);
end


%% Other Periods
for i1= 1:tr-tb
    t= tr-tb-i1+1;
    display(t)
    for i3=1:ncash
        if (i3==1)
            minc = gcash(i3,1)/5;
            maxc = 0.999*gcash(i3,1);
        else
            minc = C(i3-1,t);
            if (i3<10)
                maxc = minc + (gcash(i3,1) - gcash(i3-1,1));
            else
                mpc = max((C(i3-1,t)-C(i3-9,t))/(gcash(i3-1,1) - gcash(i3-9,1)),0.1);
                maxc = minc + mpc*(gcash(i3,1) - gcash(i3-1,1));
            end
        end
        stepc = (maxc-minc)/(nc-1);
        for i9=1:nc
            gc(i9,1)=minc+(i9-1.0)*stepc;
        end
        for i4=1:nc
            u=(1.0-delta)*(gc(i4,1)^psi_1);
            sav = gcash(i3,1)-gc(i4,1);
            for i5=1:na
                auxVV=0.0;
                for i6=1:n
                    for i8=1:n
                        for i7=1:n
                            cash_1 = riskret(i5,i8)*sav/gyp(i6,i8,t)+yh(i7,i8);
                            cash_1 = max(min(cash_1,gcash(ncash,1)),gcash(1,1));
                            int_V  = f_sc_splint(gcash(:,1),V(:,t+1),secd(:,1),ncash,cash_1);
                            auxVV  = auxVV+nweig1(i6,i7,i8)*survprob(t,1)*((int_V*gyp(i6,i8,t))^(1.0-rho));
                        end
                    end
                end
                auxV(i5,i4) = (u+delta*(auxVV^(1.0/theta)))^psi_2;
            end
        end
        vec_V = reshape(auxV,[na*nc,1]);
        [V(i3,t),pt] = max(vec_V(:,1));
        aux2 = floor((pt(1)-1)/na);
        C(i3,t) = gc(aux2+1,1);
        A(i3,t) = ga(pt(1)-na*aux2,1);
    end

    secd(:,1) = f_spline(gcash(:,1),V(:,t),ncash,1.0);

end


%% Simulations

for i1=1:floor(nsim/2)
    eps_y(1,1) = f_randn(1);
    simPY(1,i1) = eps_y(1,1)*smav;
    simPY(1,floor(nsim/2+i1)) = -eps_y(1,1)*smav;
    simGPY(1,i1) = 1.0;
    simGPY(1,floor(nsim/2+i1)) = 1.0;
    simTY(1,1) = f_randn(1);
    simY(1,i1) = exp(simTY(1,1)*smay);
    simY(1,floor(nsim/2+i1)) = exp(-simTY(1,1)*smay);

    for i2=2:tr-tb
        w = i2+tb-1;
        eps_y(1,1) = f_randn(1);
        simPY(i2,i1) = eps_y(1,1)*smav+simPY(i2-1,i1);
        simPY(i2,nsim/2+i1) = -eps_y(1,1)*smav+simPY(i2-1,i1);
        simGPY(i2,i1) = exp(gy(i2-1,1))*exp(simPY(i2,i1))/exp(simPY(i2-1,i1));
        simGPY(i2,nsim/2+i1) = exp(gy(i2-1,1))*exp(simPY(i2,nsim/2+i1))/exp(simPY(i2-1,nsim/2+i1));
        simTY(1,1) = f_randn(1);
        simY(i2,i1) = exp(simTY(1,1)*smay);
        simY(i2,nsim/2+i1) = exp(-simTY(1,1)*smay);
    end
end

for t=tr-tb+1:tn
    simY(t,:) = ret_fac;
    simGPY(t,:) = 1.0;
end

for t=1:tn
    for i1=1:floor(nsim/2)
        eps_r(1,1) = f_randn(1);
        simR(t,i1) = mu + r + sigr*eps_r(1,1);
        simR(t,nsim/2+i1) = mu + r - sigr*eps_r(1,1);
    end
end


simW(:,:) = 0.0;
for t=1:tn
    for i1=1:nsim
        simW_Y(t,i1) = simW(t,i1)/simY(t,i1);
        cash = simW(t,i1)+simY(t,i1);
        i_net_cash2 = f_ntoil(log(cash),lgcash(:,1),ncash);
        ic1 = i_net_cash2;
        ic2 = i_net_cash2+1;
        ttc = (cash-gcash(ic1,1))/(gcash(ic2,1)-gcash(ic1,1));
        ttc = max(min(1.0,ttc),0.0);
        simC(t,i1) = (1-ttc)*C(ic1,t)+ttc*C(ic2,t);
        simA(t,i1) = (1-ttc)*A(ic1,t)+ttc*A(ic2,t);
        simC(t,i1) = min(simC(t,i1),0.9999*cash);
        sav = cash-simC(t,i1);
        simS(t,i1) = simA(t,i1)*sav;
        simS(t,i1) = min(simS(t,i1),sav);
        simB(t,i1) = sav-simS(t,i1);
        if (t<tn)
            simW(t+1,i1) = (simB(t,i1)*r+simS(t,i1)*simR(t,i1))/simGPY(t+1,i1);
        end
    end
end

meanC(:,1) = simC(:,:)*ones_nsim_1(:,1)/nsim;
meanY(:,1) = simY(:,:)*ones_nsim_1(:,1)/nsim;
meanW(:,1) = simW(:,:)*ones_nsim_1(:,1)/nsim;
meanS(:,1) = simS(:,:)*ones_nsim_1(:,1)/nsim;
meanB(:,1) = simB(:,:)*ones_nsim_1(:,1)/nsim;
meanWY(:,1) = simW_Y(:,:)*ones_nsim_1(:,1)/nsim;
meanalpha(:,1) = simA(:,:)*ones_nsim_1(:,1)/nsim;
meanGPY(:,1) = simGPY(:,:)*ones_nsim_1(:,1)/nsim;

ages = (1:tn) + tb - 1;
fig = figure('Name', 'Gomes(2020)', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(ages, meanW, 'b-', 'LineWidth', 2, 'DisplayName', 'Wealth ($\bar{W}$)');
hold on;
plot(ages, meanC, 'r--', 'LineWidth', 2, 'DisplayName', 'Consumption ($\bar{C}$)');
plot(ages, meanY, 'k:', 'LineWidth', 2, 'DisplayName', 'Income ($\bar{Y}$)');
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;
title('Mean Lifecycle Profiles'); xlabel('Age'); ylabel('Value (absolute)');
grid on; xlim([tb, td]);

subplot(2,2,2);
plot(ages, meanalpha, 'm-', 'LineWidth', 2);
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5);
title('Mean Risky Asset Allocation'); xlabel('Age'); ylabel('Risky Share ($\bar{\alpha}$)', 'Interpreter', 'latex');
ylim([-0.05, 1.05]); grid on; xlim([tb, td]);

subplot(2,2,4);
plot(ages, meanWY, 'c-', 'LineWidth', 2);
xline(tr, 'k:', 'Retirement', 'LineWidth', 1.5);
title('Mean Wealth-to-Income Ratio'); xlabel('Age'); ylabel('$\bar{W} / \bar{Y}$', 'Interpreter', 'latex');
grid on; xlim([tb, td]);

sgtitle('Gomes(2020)', 'FontSize', 16, 'FontWeight', 'bold');

output_dir = 'fig';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
output_filename = fullfile(output_dir, 'gomes.png');
print(fig, output_filename, '-dpng', '-r300');