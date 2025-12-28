%% Life-Cycle Consumption/Portfolio Choice Problem

clear all;
close all;
warning off
%% Variable Definitions
filename = string(zeros(80,1));

tb     = 20; % initial starting age
tr     = 66; % retire age
td     = 100; % death age
tn     = td-tb+1; % total period

% na     = 11; % asset grid num
% ncash  = 11; % cash grid num
% n      = 5; % five stochastic state
% nc     = 101; % consumption grid
na     = 21;
ncash  = 21;
n      = 3;
nc     = 51;


maxcash     = 200.0;
mincash     = 0.25;
% 基础收入f的系数
aa          = -2.170042+2.700381;
b1          = 0.16818;
b2          = -0.0323371/10;
b3          = 0.0019704/100;

% ret_fac     = 0; % 固定养老金
ret_fac     = 0.68212; % 固定养老金
smay        = sqrt(0.0738); % 0.1; % 白噪声shock的标准差
smav        =  sqrt(0.016);%0.1; % 持续shock的标准差
corr_v      = 0.0; % 工资收入白噪声与风险收益随机项的相关性
corr_y      = 0.0; % 工资收入AR(1)随机项与风险收益随机项的相关性
rho         = 10; %3.84; % Epstein-Zin的相对风险规避系数
delta       = 0.96; %;0.97; % Epstein-Zin的discount factor, 对应Gomes(2020)中的beta
% psi         = 0.3; % Epstein-Zin的跨期替代弹性, 对应Gomes(2020)中的psi
rf           = 0.02; % 1.015; % 无风险总收入
mu          = 0.06; % 0.04; % 超额收益
sigr        = 0.157; % 0.2; % 风险资产收益率的标准差

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
c           = zeros(ncash,tn);
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
% rf = mu + Rf + epsilon*sigma, epsilon ~ N(0,1)

% realized value of epsilon
% grid(1,1) = -2.85697001387280;
% grid(2,1) = -1.35562617997427;
% grid(3,1) = 0.00000000000000;
% grid(4,1) = 1.35562617997426;
% grid(5,1) = 2.85697001387280;
%
% % probability
% weig(1,1) = 0.01125741132772;
% weig(2,1) = 0.22207592200561;
% weig(3,1) = 0.53333333333333;
% weig(4,1) = 0.22207592200561;
% weig(5,1) = 0.01125741132772;
gamma00 = 0; % AR自相关系数
mew = 0; % AR的常数项
sigma = 1; % 白噪声的标准差
tauchenoptions.parallel=0;
[grid,weig] = discretizeAR1_Tauchen(mew,gamma00,sigma,n,2,tauchenoptions);
farmertodaoptions.parallel=0;
% [grid,weig] = discretizeAR1_FarmerToda(mew,gamma00,sigma,n,farmertodaoptions);
weig = diag(weig);

%% Conditional Survival Probabilities

% survprob(1,1)  = 0.99845;
% survprob(2,1)  = 0.99839;
% survprob(3,1)  = 0.99833;
% survprob(4,1)  = 0.9983;
% survprob(5,1)  = 0.99827;
% survprob(6,1)  = 0.99826;
% survprob(7,1)  = 0.99824;
% survprob(8,1)  = 0.9982;
% survprob(9,1)  = 0.99813;
% survprob(10,1) = 0.99804;
% survprob(11,1) = 0.99795;
% survprob(12,1) = 0.99785;
% survprob(13,1) = 0.99776;
% survprob(14,1) = 0.99766;
% survprob(15,1) = 0.99755;
% survprob(16,1) = 0.99743;
% survprob(17,1) = 0.9973;
% survprob(18,1) = 0.99718;
% survprob(19,1) = 0.99707;
% survprob(20,1) = 0.99696;
% survprob(21,1) = 0.99685;
% survprob(22,1) = 0.99672;
% survprob(23,1) = 0.99656;
% survprob(24,1) = 0.99635;
% survprob(25,1) = 0.9961;
% survprob(26,1) = 0.99579;
% survprob(27,1) = 0.99543;
% survprob(28,1) = 0.99504;
% survprob(29,1) = 0.99463;
% survprob(30,1) = 0.9942;
% survprob(31,1) = 0.9937;
% survprob(32,1) = 0.99311;
% survprob(33,1) = 0.99245;
% survprob(34,1) = 0.99172;
% survprob(35,1) = 0.99091;
% survprob(36,1) = 0.99005;
% survprob(37,1) = 0.98911;
% survprob(38,1) = 0.98803;
% survprob(39,1) = 0.9868;
% survprob(40,1) = 0.98545;
% survprob(41,1) = 0.98409;
% survprob(42,1) = 0.9827;
% survprob(43,1) = 0.98123;
% survprob(44,1) = 0.97961;
% survprob(45,1) = 0.97786;
% survprob(46,1) = 0.97603;
% survprob(47,1) = 0.97414;
% survprob(48,1) = 0.97207;
% survprob(49,1) = 0.9697;
% survprob(50,1) = 0.96699;
% survprob(51,1) = 0.96393;
% survprob(52,1) = 0.96055;
% survprob(53,1) = 0.9569;
% survprob(54,1) = 0.9531;
% survprob(55,1) = 0.94921;
% survprob(56,1) = 0.94508;
% survprob(57,1) = 0.94057;
% survprob(58,1) = 0.9357;
% survprob(59,1) = 0.93031;
% survprob(60,1) = 0.92424;
% survprob(61,1) = 0.91717;
% survprob(62,1) = 0.90922;
% survprob(63,1) = 0.90089;
% survprob(64,1) = 0.89282;
% survprob(65,1) = 0.88503;
% survprob(66,1) = 0.87622;
% survprob(67,1) = 0.86576;
% survprob(68,1) = 0.8544;
% survprob(69,1) = 0.8423;
% survprob(70,1) = 0.82942;
% survprob(71,1) = 0.8154;
% survprob(72,1) = 0.80002;
% survprob(73,1) = 0.78404;
% survprob(74,1) = 0.76842;
% survprob(75,1) = 0.75382;
% survprob(76,1) = 0.73996;
% survprob(77,1) = 0.72464;
% survprob(78,1) = 0.71057;
% survprob(79,1) = 0.6961;
% survprob(80,1) = 0.6809;
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
%% Output Files - Names

% for t=1:tn-1
%     if t<10
%         filename(t,1) = string('year0')+string(t);
%     else
%         filename(t,1) = string('year')+string(t);
%     end
% end
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

% theta = (1.0-rho)/(1.0-1.0/psi);
% psi_1 = 1.0-1.0/psi;
% psi_2 = 1.0/psi_1;
%% Grids for the State Variables and for Portfolio Rule

% grids for risky asset holdings
for i1=1:na
    ga(i1,1)=(na-i1)/(na-1.0);
end

% 每个state和holdings下，投资组合(总)收益率
for i5=1:na
    for i8=1:n
        riskret(i5,i8)=rf*(1-ga(i5,1))+gret(i8,1)*ga(i5,1);
    end
end

% cash-on-hand的grid
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
    grid2(:,1) = grid(i1,1)*corr_y+grid(:,1).*ones_n_1(:,1)*(1-corr_y^2)^(0.5);
    yh(1:n,i1) = exp(grid2(:,1)*smay); % 白噪声随机项u_t的grid
end

for i1=1:n
    grid2(:,1) = grid(i1,1)*corr_v+grid(:,1).*ones_n_1(:,1)*(1-corr_v^2)^(0.5);
    yp(:,i1) = grid2(:,1)*smav;  % permanent labor p_t 的随机游走项白噪声z_t的grid
end

% 随年龄变化的基础收入 随年龄先递增后递减
for i1=tb:tr
    f_y(i1-tb+1,1) = exp(aa+b1*i1+b2*i1^2+b3*i1^3); % f(t,theta)
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
C(:,tn) = 1; % 最后一期的现金全部用来消费
A(:,tn) = 0.0; % 最后一期的投资全部为零

for i1=1:ncash
    V(i1,tn) = (C(i1,tn)*gcash(i1))^(1-rho)/(1-rho); % 最后一期的值函数（Epstein-Zin效用函数）
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
        x0 = [0.1,0.1];
        lb = [0,0];
        ub = [1,1];
        options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
        policy = fmincon(@(x) cocco_fun_valuefunc_retired(x,gcash(i3),V(:,t+1),...
            gret,rf,ret_fac,rho,delta,gcash,survprob(t,1),weig),...
            x0,[],[],[],[],lb,ub,[],options);
        C(i3,t) = policy(1);
        A(i3,t) = policy(2);
        V(i3,t) = -cocco_fun_valuefunc_retired(policy,gcash(i3),V(:,t+1),...
            gret,rf,ret_fac,rho,delta,gcash,survprob(t,1),weig);
        % maxc = 1;
        % minc = 0;
        % % 确定消费比例的grid
        % stepc=(maxc-minc)/(nc-1);
        % for i9=1:nc
        %     gc(i9,1)=minc+(i9-1.0)*stepc;
        % end
        % % 开始迭代
        % % ========================  *控制变量c ===================================
        % for i4=1:nc
        %     u=(1.0-delta)*((gc(i4,1)*gcash(i3,1))^psi_1) ;
        %     sav = gcash(i3,1)*(1-gc(i4,1));
        %     % ========================  *控制变量a ===================================
        %     for i5=1:na % 资产权重grid
        %         auxVV=0.0;
        %         %           for i6=1:n 无工资的不确定性
        %         %              for i7=1:n
        %         % ========================  *外生状态变量:只有资产收益的随机项 =============
        %         for i8=1:n % 随机资产收益率的grid
        %             cash_1 = riskret(i5,i8)*sav + ret_fac; % 下期的cash-on-hand
        %             cash_1 = max(min(cash_1,gcash(ncash,1)),gcash(1,1)); % 确保下期cash在gcash的[min,max]范围内
        %             % log_cash(i3,:) = cash_1;
        %             % int_V  = f_sc_splint(gcash(:,1),V(:,t+1),secd(:,1),ncash,cash_1); % cubic spline 插值计算V(下期cash)
        %             int_V = interp1(gcash(:,1), V(:,t+1), cash_1, 'cubic'); % matlab自带的插值，无需利用f_spline.m手动计算二次导数
        %             auxVV= auxVV + weig(i8,1)*survprob(t,1)*(int_V^(1.0-rho)); % 累加auxVV得到EV
        %         end
        %         %              end
        %         %           end
        %         % 得到所有可能的control var的grid下的V
        %         auxV(i5,i4) = (u+delta*(auxVV^(1.0/theta)))^psi_2;    %1/(1-1/psi), Epstein-Zin效用形式
        %     end
        % end
        % vec_V = reshape(auxV,[na*nc,1]); % 寻找最大的V，作为t期的V
        % [V(i3,t),pt] = max(vec_V(:,1));
        % aux2 = floor((pt(1)-1)/na);
        % % 最大V对应控制变量grid的值（policy function）
        % C(i3,t) = gc(aux2+1,1)*gcash(i3,1); % 消费实际值
        % c(i3,t) = gc(aux2+1,1); % 消费比例
        % A(i3,t) = ga(pt(1)-na*aux2,1);
    end
end

%%
% Retirement Periods - Save to txt

% for i1 = 1:td-tr
%    t = tn-i1;
%     fileID = fopen(filename(t,1)+'.txt','w');
%     for i5=1:ncash
%         fprintf(fileID,'%12.8f %12.8f\rf\n',A(i5,t),gcash(i5,1));
%     end
%     for i5=1:ncash
%         fprintf(fileID,'%12.8f %12.8f\rf\n',C(i5,t),gcash(i5,1));
%     end
%     for i5=1:ncash
%         fprintf(fileID,'%12.8f %12.8f\rf\n',V(i5,t),gcash(i5,1));
%     end
%     fclose(fileID);
% end
%% Other Periods（退休前）

% ======================== *内生状态变量t =================================
for i1= 1:tr-tb
    t= tr-tb-i1+1;
    display(t)
    % ======================== *内生状态变量cash ==============================
    for i3=1:ncash
        x0 = [0.1,0.1];
        lb = [0,0];
        ub = [1,1];
        options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
        policy = fmincon(@(x) cocco_fun_valuefunc_work(x,gcash(i3),gyp(:,:,t),V(:,t+1),...
            yh,gret,rf,rho,delta,gcash,survprob(t,1),n,nweig1),...
            x0,[],[],[],[],lb,ub,[],options);
        C(i3,t) = policy(1);
        A(i3,t) = policy(2);
        V(i3,t) = -cocco_fun_valuefunc_work(policy,gcash(i3),gyp(:,:,t),V(:,t+1),...
            yh,gret,rf,rho,delta,gcash,survprob(t,1),n,nweig1);
        % maxc = 1;
        % minc = 0;
        % stepc = (maxc-minc)/(nc-1);
        % for i9=1:nc
        %     gc(i9,1)=minc+(i9-1.0)*stepc;
        % end
        % % 开始迭代
        % % ========================  *控制变量c ===================================
        % for i4=1:nc
        %     u=(1.0-delta)*((gc(i4,1)*gcash(i3,1))^psi_1);
        %     sav = gcash(i3,1)*(1-gc(i4,1));
        %     % ========================  *控制变量a ==================================
        %     for i5=1:na
        %         auxVV=0.0;
        %         % ========================  *外生状态变量:3个随机项 =====================
        %         for i6=1:n
        %             for i8=1:n
        %                 for i7=1:n
        %                     cash_1 = riskret(i5,i8)*sav/gyp(i6,i8,t) + yh(i7,i8); % 这里的cash是scaled后的cash
        %                     cash_1 = max(min(cash_1,gcash(ncash,1)),gcash(1,1));
        %                     int_V  = f_sc_splint(gcash(:,1),V(:,t+1),secd(:,1),ncash,cash_1);
        %                     % int_V = interp1(gcash(:,1), V(:,t+1), cash_1, 'cubic'); % matlab自带的插值，无需利用f_spline.m手动计算二次导数
        %                     auxVV  = auxVV+nweig1(i6,i7,i8)*survprob(t,1)*((int_V*gyp(i6,i8,t))^(1.0-rho)); % 新的scaled value function
        %                 end
        %             end
        %         end
        %         auxV(i5,i4) = (u+delta*(auxVV^(1.0/theta)))^psi_2;
        %     end
        % end
        % vec_V = reshape(auxV,[na*nc,1]);
        % [V(i3,t),pt] = max(vec_V(:,1));
        % aux2 = floor((pt(1)-1)/na);
        % C(i3,t) = gc(aux2+1,1)*gcash(i3,1);  % 消费实际值
        % c(i3,t) = gc(aux2+1,1); % 消费比例
        % A(i3,t) = ga(pt(1)-na*aux2,1);
    end
end

% currentDir = fileparts(mfilename('fullpath')); % 获取当前m文件所在目录
% % save([currentDir,'\result\Cocco_', strrep(datestr(datetime('now')), ':', '-'),'.mat'])
% save([currentDir,'\result\Cocco.mat'])

%% 数值模拟
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
simW(:,:) = 1;
for t=1:tn
    for i1=1:nsim
        simW_Y(t,i1) = simW(t,i1)/simY(t,i1); % 上期财富-本期工资收入比
        cash(t,i1) = simW(t,i1)+simY(t,i1); % cash-on-hand
        i_net_cash2 = f_ntoil_modified(log(cash(t,i1)),lgcash(:,1),ncash);  % 寻找value在cash grid中最接近的位置，超出的话就返回边界的位置
        ic1 = i_net_cash2;
        ic2 = i_net_cash2+1;
        ttc = (cash(t,i1)-gcash(ic1,1))/(gcash(ic2,1)-gcash(ic1,1));
        ttc = max(min(1.0,ttc),0.0); % 确保插值权重在0和1之间
        simC(t,i1) = (1-ttc)*C(ic1,t)+ttc*C(ic2,t); % 给定cash，在最优的消费比例之间插值
        simc(t,i1) = cash(t,i1)*simC(t,i1);
        % simC(t,i1) = interp2(gcash,C(:,t),gcash,'spline')';
        simA(t,i1) = (1-ttc)*A(ic1,t)+ttc*A(ic2,t); % 给定cash，在最优的权重之间插值
        % simA(t,i1) = interp2(gcash,A(:,t),gcash,'spline')';


        sav = cash(t,i1)*(1-simC(t,i1)); % 用于投资的金额
        simS(t,i1) = simA(t,i1)*sav; % 风险投资额
        simS(t,i1) = min(simS(t,i1),sav);
        simB(t,i1) = sav-simS(t,i1); % 无风险投资额
        if (t<tn)
            simW(t+1,i1) = (simB(t,i1)*rf+simS(t,i1)*simR(t,i1))/simGPY(t+1,i1); % 下期一开始的(scaled)财富（拿到下期工资收入前）
        end
    end
end

% 多次模拟path下变量平均值
meanC(:,1) = simC(:,:)*ones_nsim_1(:,1)/nsim;
meanc(:,1) = simc(:,:)*ones_nsim_1(:,1)/nsim;
meanY(:,1) = simY(:,:)*ones_nsim_1(:,1)/nsim;
meanW(:,1) = simW(:,:)*ones_nsim_1(:,1)/nsim;
meanS(:,1) = simS(:,:)*ones_nsim_1(:,1)/nsim;
meanB(:,1) = simB(:,:)*ones_nsim_1(:,1)/nsim;
meanWY(:,1) = simW_Y(:,:)*ones_nsim_1(:,1)/nsim;
meanalpha(:,1) = simA(:,:)*ones_nsim_1(:,1)/nsim;
meanGPY(:,1) = simGPY(:,:)*ones_nsim_1(:,1)/nsim;

plot(meanc);
hold on
plot(meanalpha);
hold on
% plot(meanW);
legend('c','\alpha')



