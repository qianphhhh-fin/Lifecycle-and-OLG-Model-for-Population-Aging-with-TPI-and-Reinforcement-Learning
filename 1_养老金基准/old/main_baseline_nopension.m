%% Life-Cycle Consumption/Portfolio Choice Problem

clear all;
close all;
warning off
%% Variable Definitions
filename = string(zeros(80,1));

tb = 18; %96; % ; % 初始开始的年纪
tr = 61; %98; % 61; % 退休年龄
td = 100; % 死亡年龄
tn = td-tb+1; % 总期数
% parfor循环设置
% numWorkers = 8;
% parpool(16);
% 控制变量grid数量

% na  = 11; % (0,1), 风险资产
% nc = 11; % R+, 消费
% nq = 11; % {0,1}, 是否购买养老保险
% 状态变量grid数量

ncash  = 51; % 手中现金
n = 5; % 外生随机冲击的grid数量
% 外生参数

% 基础收入f的系数
aa          = (-2.170042+2.700381);
b1          = 0.16818;
b2          = -0.0323371/10;
b3          = 0.0019704/100;
% aa = -1.40432388e+01;
% b1 = 1.42492731e+00;
% b2 = -4.55146961e-02;
% b3 = 6.08573382e-04;
% b4 = -2.88233535e-06;
% aa          = 8.5473;
% b1          = 0.0809;
% b2          = -0.0004;
% b3          = -0.00001113;


% 养老金相关
ret_fac     = 0.6827; % 退休后固定支付的工资（养老金）
% ret_fac     = 3; % 固定养老金
% kappa = 0.1; % 将养老金固定投资于风险资产的比例（不由家庭决定）
% zeta = 1/(td-tr+1); % 在退休期，每年取出养老基金的zeta支付给消费者
pension_pct = 0; % 工资扣除缴纳养老保险的比例

smay = sqrt(0.169993); % ; % 白噪声shock的标准差
smav = sqrt(0.112572); % 持续shock的标准差
corr_z_epsilon = 0.0; % 工资收入白噪声与风险收益随机项的相关性
corr_u_epsilon = 0.0; % 工资收入AR(1)随机项与风险收益随机项的相关性

gamma = 3.84; % Epstein-Zin的相对风险规避系数
beta = 0.95; % Epstein-Zin的discount factor
psi = 0.15; % Epstein-Zin的跨期替代弹性, 越高越看重当期消费
rf = 1.02; % 无风险总收入
mu = 0.04; % 超额收益
sigr = 0.27; % 风险资产收益率的标准差
% r_pension = 1.04; % 养老基金的收益率

survprob    = zeros(tn-1,1);
beta2      = zeros(tn-1,1);
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
%% Approximation to Normal Distribution 模拟一个服从正态分布N(0,1)的白噪声.离散化
% R = mu + Rf + epsilon*sigma, epsilon ~ N(0,1)
gamma00 = 0; % AR自相关系数
mew = 0; % AR的常数项
sigma = 1; % 白噪声的标准差
tauchenoptions.parallel=0;
[grid,weig] = discretizeAR1_Tauchen(mew,gamma00,sigma,n,2,tauchenoptions);
farmertodaoptions.parallel=0;
% [grid,weig] = discretizeAR1_FarmerToda(mew,gamma00,sigma,n,farmertodaoptions);
weig = diag(weig);

%% Conditional Survival Probabilities

survprob(1,1)=0.99975;
survprob(2,1)=0.99974;
survprob(3,1)=0.99973;
survprob(4,1)=0.99972;
survprob(5,1)=0.99971;
survprob(6,1)=0.99969;
survprob(7,1)=0.99968;
survprob(8,1)=0.99966;
survprob(9,1)=0.99963;
survprob(10,1)=0.99961;
survprob(11,1)=0.99958;
survprob(12,1)=0.99954;
survprob(13,1)=0.99951;
survprob(14,1)=0.99947;
survprob(15,1)=0.99943;
survprob(16,1)=0.99938;
survprob(17,1)=0.99934;
survprob(18,1)=0.99928;
survprob(19,1)=0.99922;
survprob(20,1)=0.99916;
survprob(21,1)=0.99908;
survprob(22,1)=0.999;
survprob(23,1)=0.99891;
survprob(24,1)=0.99881;
survprob(25,1)=0.9987;
survprob(26,1)=0.99857;
survprob(27,1)=0.99843;
survprob(28,1)=0.99828;
survprob(29,1)=0.99812;
survprob(30,1)=0.99794;
survprob(31,1)=0.99776;
survprob(32,1)=0.99756;
survprob(33,1)=0.99735;
survprob(34,1)=0.99713;
survprob(35,1)=0.9969;
survprob(36,1)=0.99666;
survprob(37,1)=0.9964;
survprob(38,1)=0.99613;
survprob(39,1)=0.99583;
survprob(40,1)=0.99551;
survprob(41,1)=0.99515;
survprob(42,1)=0.99476;
survprob(43,1)=0.99432;
survprob(44,1)=0.99383;
survprob(45,1)=0.9933;
survprob(46,1)=0.9927;
survprob(47,1)=0.99205;
survprob(48,1)=0.99133;
survprob(49,1)=0.99053;
survprob(50,1)=0.98961;
survprob(51,1)=0.98852;
survprob(52,1)=0.98718;
survprob(53,1)=0.98553;
survprob(54,1)=0.98346;
survprob(55,1)=0.98089;
survprob(56,1)=0.97772;
survprob(57,1)=0.97391;
survprob(58,1)=0.96943;
survprob(59,1)=0.96429;
survprob(60,1)=0.95854;
survprob(61,1)=0.95221;
survprob(62,1)=0.94537;
survprob(63,1)=0.93805;
survprob(64,1)=0.93027;
survprob(65,1)=0.92202;
survprob(66,1)=0.91327;
survprob(67,1)=0.90393;
survprob(68,1)=0.89389;
survprob(69,1)=0.88304;
survprob(70,1)=0.87126;
survprob(71,1)=0.85846;
survprob(72,1)=0.84452;
survprob(73,1)=0.82935;
survprob(74,1)=0.81282;
survprob(75,1)=0.79485;
survprob(76,1)=0.77543;
survprob(77,1)=0.75458;
survprob(78,1)=0.7324;
survprob(79,1)=0.70893;
survprob(80,1)=0.68424;
survprob(81,1)=0.68424;
survprob(82,1)=0.68424;



% survprob = ones(80,1); % 无死亡概率

%% Additional Computations
% risky (total) return
for i1=1:n
    gret(i1,1) = rf+mu+grid(i1,1)*sigr;
end

% 外生随机性的概率
for i6=1:n
    for i7=1:n
        for i8=1:n
            nweig1(i6,i7,i8) = weig(i6,1)*weig(i7,1)*weig(i8,1);
        end
    end
end

theta = (1.0-gamma)/(1.0-1.0/psi);
psi_1 = 1.0-1.0/psi;
psi_2 = 1.0/psi_1;
%% Grids for the State Variables and for Portfolio Rule


% cash-on-hand的grid
maxcash = 100;
mincash = 0.25;
l_maxcash = log(maxcash);
l_mincash = log(mincash);
stepcash = (l_maxcash-l_mincash)/(ncash-1);

for i1=1:ncash
    lgcash(i1,1)=l_mincash+(i1-1.0)*stepcash;
end
for i1=1:ncash
    gcash(i1,1)=exp(lgcash(i1,1));
end
%% Labor Income (normalization)

for i1=1:n
    grid2(:,1) = grid(i1,1)*corr_u_epsilon+grid(:,1).*ones_n_1(:,1)*(1-corr_u_epsilon^2)^(0.5);
    yh(1:n,i1) = exp(grid2(:,1)*smay); % 白噪声随机项u_t的grid
end

for i1=1:n
    grid2(:,1) = grid(i1,1)*corr_z_epsilon+grid(:,1).*ones_n_1(:,1)*(1-corr_z_epsilon^2)^(0.5);
    yp(:,i1) = grid2(:,1)*smav;  % permanent labor p_t 的随机游走项白噪声z_t的grid
end

% 随年龄变化的基础收入 随年龄先递增后递减
for i1=tb:tr
    f_y(i1-tb+1,1) = exp((aa+b1*i1+b2*i1^2+b3*i1^3)); % f(t,theta)+b4*i1^4
end

% working age
for i1=tb:tr-1
    gy(i1-tb+1,1) = f_y(i1-tb+2,1)/f_y(i1-tb+1,1)-1.0; % 对基础工资收入进行normalization
    %     gy(i1-tb+1,1) = f_y(i1-tb+2,1)-f_y(i1-tb+1,1); % 将整个优化问题scale后剩下要除的term
    for i2=1:n
        gyp(:,i2,i1-tb+1) = exp(gy(i1-tb+1,1)*ones_n_1(:,1)+yp(:,i2)); % normalized工资收入+持久收入shock
    end
end

% retirement age
for i1=tr-tb+1:tn-1
    for i2=1:n
        gyp(:,i2,i1) = exp(0.0*ones_n_1(:,1));
    end
end
%% Terminal Period

% policy function
C(:,tn) = gcash(:,1); % 最后一期的现金全部用来消费

A(:,tn) = 0.0; % 最后一期的投资全部为零
delta =0.2;
for i1=1:ncash
    % V(i1,tn) = C(i1,tn)^(1-gamma)/(1-gamma); % 最后一期的值函数（Epstein-Zin效用函数）
    % V(i1,tn) = (C(i1,tn))^(1-gamma)/(1-gamma);
    % if C(i1,tn)<delta
    %     c = -gamma/2*delta^(-gamma-1);
    %     b = delta^(-gamma)-2*c*delta;
    %     a = delta^(1-gamma)/(1-gamma)-b*delta-c*delta^2;
    %     V(i1,tn)= a+b*C(i1,tn)+c*C(i1,tn)^2;
    % else
        V(i1,tn) = -(C(i1,tn))^(1-gamma);
    % V(i1,tn) = log(C(i1,tn));
    % % end
    % V(i1,tn) = (C(i1,tn))^(1-gamma)/(1-gamma);
    % V(i1,tn) = log(C(i1,tn));
    % if V(i1,tn)<-2
    %     V(i1,tn)=-2;
    % end
    % V(i1,tn) = log(1+C(i1,tn)^(1-gamma)/(1-gamma));
end

% ================ VFI ============================
%% Retirement Periods

% ======================== *内生状态变量t =================================
for i1 = 1:td-tr
    t = tn-i1;
    display(t)
    %    minmaxcash = nan(ncash,2);
    % ======================== *内生状态变量cash ==============================
    for i3=1:ncash

        x0 = [0.2,0.2];
        lb = [0,0];
        ub = [0.999*gcash(i3,1),1];
        if i3~=1
            lb = [C(i3-1,t),0];
            % ub = [0.999*gcash(i3,1),A(i3-1,t)];
        end
        options = optimoptions('fmincon','Display','off','Algorithm','sqp');
        policy = fmincon(@(x) cocco_fun_valuefunc_retired(x,gcash(i3),V(:,t+1),...
            gret,rf,ret_fac,gamma,beta,psi_1,psi_2,theta,gcash,survprob(t,1),weig),...
            x0,[],[],[],[],lb,ub,[],options);
        C(i3,t) = policy(1);
        A(i3,t) = policy(2);
        V(i3,t) = -cocco_fun_valuefunc_retired(policy,gcash(i3),V(:,t+1),...
            gret,rf,ret_fac,gamma,beta,psi_1,psi_2,theta,gcash,survprob(t,1),weig);
        % -cocco_fun_valuefunc_retired([C(i3,t),A(i3,t)],gcash(i3),V(:,t+1),...
        %     gret,r,ret_fac,gamma,beta,psi_1,psi_2,theta,gcash,survprob(t,1),weig);
    end
end


%% Other Periods（退休前）

% ======================== *内生状态变量t =================================
for i1= 1:tr-tb
    t= tr-tb-i1+1;
    display(t)
    % ======================== *内生状态变量cash ==============================
    for i3=1:ncash

        x0 = [0.1,0.1];
        lb = [0,0];
        ub = [0.999*gcash(i3,1),1];
        if i3~=1
            lb = [C(i3-1,t),0];
            % ub = [0.999*gcash(i3,1),A(i3-1,t)];
        end
        options = optimoptions('fmincon','Display','off','Algorithm','sqp');
        policy = fmincon(@(x) cocco_fun_valuefunc_work(x,gcash(i3),gyp(:,:,t),V(:,t+1),...
            yh,gret,rf,gamma,beta,psi_1,psi_2,theta,gcash,survprob(t,1),n,nweig1),...
            x0,[],[],[],[],lb,ub,[],options);
        C(i3,t) = policy(1);
        A(i3,t) = policy(2);
        V(i3,t) = -cocco_fun_valuefunc_work(policy,gcash(i3),gyp(:,:,t),V(:,t+1),...
            yh,gret,rf,gamma,beta,psi_1,psi_2,theta,gcash,survprob(t,1),n,nweig1);


    end


end

% % 
writematrix(C, 'result_cocco_matlab\C.xlsx','WriteMode','overwrite');  % 基础数值数据
writematrix(A, 'result_cocco_matlab\A.xlsx','WriteMode','overwrite');  % 基础数值数据
writematrix(V, 'result_cocco_matlab\V.xlsx','WriteMode','overwrite');  % 基础数值数据
writematrix(gcash, 'result_cocco_matlab\gcash.xlsx','WriteMode','overwrite');  % 基础数值数据
%% 数值模拟
currentDir = fileparts(mfilename('fullpath')); % 获取当前m文件所在目录
% save([currentDir,'\result\Cocco_', strrep(datestr(datetime('now')), ':', '-'),'.mat'])
% save([currentDir,'\result\baseline_nopension.mat'])


% 1、模拟生成labor income

for i1=1:floor(nsim/2) % 另外一半模拟完全对称

    % working period第一期
    eps_y(1,1) = f_randn(1);% N(0,1)
    simPY(1,i1) = eps_y(1,1)*smav; % 初始的p
    simPY(1,floor(nsim/2+i1)) = -eps_y(1,1)*smav;
    simGPY(1,i1) = 1.0;
    simGPY(1,floor(nsim/2+i1)) = 1.0;
    simTY(1,1) = f_randn(1);
    simY(1,i1) = exp(simTY(1,1)*smay);
    simY(1,floor(nsim/2+i1)) = exp(-simTY(1,1)*smay);

    % working period第2期~退休
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

% 退休期
for t=tr-tb+1:tn
    simY(t,:) = ret_fac;
    simGPY(t,:) = 1.0;
end

% 2、模拟风险投资的收益率
for t=1:tn
    for i1=1:floor(nsim/2)
        eps_r(1,1) = f_randn(1);
        simR(t,i1) = mu + rf + sigr*eps_r(1,1);
        simR(t,nsim/2+i1) = mu + rf - sigr*eps_r(1,1);
    end
end

% 从第一期开始迭代，得到各控制变量的值
simW(:,:) = 0;
for t=1:tn
    for i1=1:nsim
        simW_Y(t,i1) = simW(t,i1)/simY(t,i1); % 上期财富-本期工资收入比
        cash(t,i1) = simW(t,i1)+simY(t,i1); % cash-on-hand
        i_net_cash2 = f_ntoil_modified(log(cash(t,i1)),lgcash(:,1),ncash);  % 寻找value在cash grid中最接近的位置，超出的话就返回边界的位置
        ic1 = i_net_cash2;
        ic2 = i_net_cash2+1;
        ttc = (cash(t,i1)-gcash(ic1,1))/(gcash(ic2,1)-gcash(ic1,1));
        ttc = max(min(1.0,ttc),0.0); % 确保插值权重在0和1之间
        simC(t,i1) = interp1(gcash,C(:,t),cash(t,i1),'spline');
        % simC(t,i1) = simC(t,i1)/cash(t,i1);
        simA(t,i1) = interp1(gcash,A(:,t),cash(t,i1),'spline');
        
        % simC(t,i1) = (1-ttc)*C(ic1,t)+ttc*C(ic2,t); % 给定cash，在最优的消费之间插值
        % simA(t,i1) = (1-ttc)*A(ic1,t)+ttc*A(ic2,t); % 给定cash，在最优的权重之间插值
        % simC(t,i1) = interp1(gcash,C(:,t),cash(t,i1),'spline')';
        %  simA(t,i1) = interp1(gcash,A(:,t),cash(t,i1),'spline')';

        simC(t,i1) = max(min(simC(t,i1),0.9999*cash(t,i1)),0);
        % reward(t,i1) = -(simC(t,i1))^(1-gamma)*50+10;
        simA(t,i1) = min(max(simA(t,i1),0),1);
        simC_pct(t,i1) = simC(t,i1)/cash(t,i1);
        sav = cash(t,i1)-simC(t,i1); % 用于投资的金额
        simS(t,i1) = simA(t,i1)*sav; % 风险投资额
        simS(t,i1) = min(simS(t,i1),sav);
        simB(t,i1) = sav-simS(t,i1); % 无风险投资额
        
        simC(t,i1) = simC(t,i1)/cash(t,i1);
        if (t<tn)
            simW(t+1,i1) = (simB(t,i1)*rf+simS(t,i1)*simR(t,i1))/simGPY(t+1,i1); % 下期一开始的(scaled)财富（拿到下期工资收入前）
        end
    end
end

% 多次模拟path下变量平均值
% meanR = reward.*[1:tn]*gamma
meanC(:,1) = mean(simC(:,:),2);
meanC_pct = mean(simC_pct(:,:),2);
meanY(:,1) = simY(:,:)*ones_nsim_1(:,1)/nsim;
meanW(:,1) = simW(:,:)*ones_nsim_1(:,1)/nsim;
meanS(:,1) = simS(:,:)*ones_nsim_1(:,1)/nsim;
meanB(:,1) = simB(:,:)*ones_nsim_1(:,1)/nsim;
meanWY(:,1) = simW_Y(:,:)*ones_nsim_1(:,1)/nsim;
meanalpha(:,1) = simA(:,:)*ones_nsim_1(:,1)/nsim;
meanGPY(:,1) = simGPY(:,:)*ones_nsim_1(:,1)/nsim;

plot(meanC_pct);
hold on
plot(meanalpha);
hold on
% plot(meanW);
legend('c','\alpha')

% save([currentDir,'\result_analysis\simu_result\baseline_nopension.mat'])


