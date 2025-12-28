% Huggett(1996)模型复现 - 生命周期经济中的财富分布
%
% 此代码求解过程比必要的时间更长，因为我希望在复现过程中对一般均衡的求解步骤进行多重检查。
%
% 我遵循Huggett(1996)的符号表示，除了时期数量(和退休年龄)，我用J(JR)表示而不是N。
%
% 一个内生变量：资产
% 一个随机外生变量：收入冲击
% 年龄

% 在服务器上运行所需的几行代码
addpath(genpath('./MatlabToolkits/'))

CheckUniquenessOfGE=0; % 如果=1，将在GE周围使用价格网格检查其他可能的GE。会大幅增加运行时间。

Params.J=79; % 年龄从20岁到98岁。

% 使用的网格大小
n_a=20;
% n_z=19; % 收入（这19个点被硬编码到z_grid和pi_z中，这样做是因为Huggett的设置方式）
N_j=Params.J; % 有限视野中的周期数

%% 声明模型参数
% 注意r, w和G将在一般均衡中确定，所以这些只是初始猜测。

% 偏好参数
Params.beta=1.011; % 个体折现未来的比率（大于1是因为"随机死亡概率"）
Params.sigma=1.5; % 风险厌恶系数

Params.bc_equalsminusw=0; % 借贷约束的两种可能值之一的指示变量，0和-w

% 人口统计学
Params.JR=46; % 65岁退休（注：实际上从未直接使用，因为它隐含在确定性收入曲线和退休福利中）
Params.n=0.012; % 人口增长率1%

% 税率
% Params.tau % 基于r确定：Params.tau=0.195/(1-Params.delta*K/Y) % 注意K/Y可以从r计算，见下文
Params.tau=0.195/(1-0.06*3); % 这里只给tau一个初始猜测
Params.theta=0.1;
% 意外遗产
% Params.T % 在GE中确定
% 退休福利：我将其设为b*bvec（在Huggett的符号中这只是b，它本身是退休状态的函数，
% 我将其分为一个标量和一个指示你是否退休的指标）。
Params.bvec=[zeros(1,Params.JR-1),ones(1,1+Params.J-Params.JR)]; % 设置为一个年龄依赖的向量，在退休年龄前为零
% 注意：bvec实际上只是退休的指标

% 生产函数
Params.A=0.895944;
Params.alpha=0.36; % Cobb-Douglas生产函数中的资本份额（我目前住在新西兰，那里的实际GDP资本份额就是这个数字 ;)
Params.delta=0.06; % 折旧率

% 生存概率，基于Jordan, C., Life contingencies, 2nd ed. (Society of Actuaries)。
% Huggett(1996)将其日期标为1975年印刷版，以下链接是1991年印刷版。但由于它仍然是第2版（首次出现于1967年）
% 似乎是正确的数字。（pdf的第342页，书的第346页，似乎是最接近的内容）
% https://vdocuments.mx/download/life-contingencies-chester-wallace-jordanpdf
Params.dj=[0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
    0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
    0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
    0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
    0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
    0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
    0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
    0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018]; % 条件死亡概率
Params.sj=1-Params.dj; % 条件生存概率。作为折现率。

% 声明年龄依赖参数。这只是创建
% 长度为J的行向量参数（VFI工具包自动
% 检测哪些参数依赖于年龄，哪些不依赖）。
% Params.ybarj=log(0.5289)+log([linspace(0.3,1.25,12),linspace(1.3,1.5,6),linspace(1.5,1.5,11),linspace(1.49,1.05,11),linspace(1,0.1,9),linspace(0.08,0,5),zeros(1,Params.J-54)]); % 确定性收入依赖于年龄
Corbae_deterministicincome = [0.0911 0.1573 0.2268 0.2752 0.3218 0.3669 0.4114 0.4559 0.4859 0.5164 0.5474 0.5786 0.6097 0.6311 0.6517 0.6711 0.6893 0.7060 0.7213 0.7355 0.7489 0.7619 0.7747 0.7783 0.7825 0.7874 0.7931 0.7994 0.7923 0.7850 0.7771 0.7679 0.7567 0.7351 0.7105 0.6822 0.6500 0.6138 0.5675 0.5183 0.4672 0.3935 0.3239 0.2596 0.1955 0.1408 0.0959 0.0604 0.0459 0.0342 0.0246 0.0165 0.0091 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % Huggett(1996)的图1显示：plot(19+(1:1:Params.J),exp(Params.ybarj))
Params.ybarj=log(Corbae_deterministicincome);
% log(0.5289)+来自Dean Corbae 2014年的作业讲义，他
% 声明这些数据来自Mark Huggett（我从朋友那里获得了Corbae的讲义文件）。
% 唯一奇怪的是乘以0.5289的标准化。我在Huggett(1996)中找不到
% 这方面的提及，但根据Corbae的作业，
% 这是为了使"工作"L与人口N的比率正确；
% 与数据中的0.5289相符。

% 随机收入，y:
Params.gamma=0.96;
Params.sigmasqepsilon=0.045;
% Params.sigmasqy=Params.sigmasqepsilon./(1-Params.gamma.^2); 
% 收入的初始分布
Params.sigmasqy1=0.38;

% Huggett使用17个状态，在+-4*sigmasqy1之间均匀分布，第18个
% 状态在6*sigmasqy1，"状态之间的转移概率是
% 通过积分正态分布下的面积计算的，条件是
% 当前状态的值。"
n_z=18;
% z_grid=[linspace(-4*sqrt(Params.sigmasqy1),4*sqrt(Params.sigmasqy1),17),6*sqrt(Params.sigmasqy1)]';
% pi_z=nan(18,18);
% % 以下几行实现转移矩阵，它们主要是从TauchenMethod()命令复制的一些代码。
% sigma=sqrt(Params.sigmasqepsilon); %e的标准差
% for ii=1:length(z_grid)
%     pi_z(ii,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2-Params.gamma*z_grid(ii),0,sigma);
%     for jj=2:(length(z_grid)-1)
%         pi_z(ii,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2-Params.gamma*z_grid(ii),0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2-Params.gamma*z_grid(ii),0,sigma);
%     end
%     pi_z(ii,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2-Params.gamma*z_grid(ii),0,sigma);
% end
% z_grid=gpuArray(z_grid);
% pi_z=gpuArray(pi_z);
% 双重检查：cumsum(pi_z,2)显示每行总和为1。

%% 一般均衡变量：给出一些初始值
GEPriceParamNames={'r','b','T'};
Params.r=0.05; % 资产利率
Params.b=1.2; % 退休人员的福利水平
Params.T=0.1; % 从意外遗产中支付的一次性转移
% 我最初设置b=0.8, T=0.6；在求解GE后切换到这些值，因为我知道它们
% 更接近真实值，这有助于加快运行速度。

% 以下内容被编码到回报函数中，以根据r的值获取w和tau的值
% KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% KdivY=(KdivL^(1-Params.alpha))/Params.A;
% Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % 工资率（每有效劳动单位）
% Params.tau=0.195/(1-Params.delta*KdivY);

% 注意G不是GEPriceParamNames的一部分，这是因为它
% 实际上只是模型的剩余部分，在实际
% 计算中不起作用。它不影响任何人的行为，所以我们不需要它来
% 解决模型，一旦我们解决了模型，如果我们想知道它是什么
% 我们可以根据政府的预算平衡方程计算它。


%% 网格
% w随r递减，所以假设r不会低于零，我们知道w不会高于：
maxw=Params.A*(1-Params.alpha)*((((0+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1)))^Params.alpha);
maxa=250;
a_grid=[linspace(-maxw,15,ceil(n_a/2))'; linspace(15,maxa,n_a-ceil(n_a/2))']; % 选择15是因为这无论如何都大于平均值
a_grid(ceil(n_a/2))=0; %有两个15值，用0替换一个
a_grid=sort(a_grid); % 双重检查：length(unique(a_grid))==n_a
% 注意借贷约束是在回报函数内部强制执行的。
% 这意味着一些网格点被浪费了，但更加整洁。

%% 现在，创建回报函数
DiscountFactorParamNames={'beta','sj'};
 
ReturnFn=@(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw) Huggett1996_ReturnFn(aprime,a,z,sigma,r,ybarj,theta,b,bvec,T,delta,alpha,A,bc_equalsminusw)
% 对于回报函数，第一个输入必须是（任何决策变量），下一期内生
% 状态，本期内生状态（任何外生冲击）。之后是任何参数。

vfoptions.verbose=0;
vfoptions.policy_forceintegertype=2; % 策略没有被视为整数（其中一个元素与整数相差10^(-15)）


%% 个体年龄分布
% 几乎所有OLG模型都包含某种人口增长，也许
% 还有一些其他因素，创造了不同年龄的加权，需要
% 用于计算稳态分布和聚合变量。
Params.mewj=ones(1,Params.J);
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj./sum(Params.mewj);

simoptions.ncores=feature('numcores'); % CPU核心数
simoptions.nsims=4*10^5; % 默认值为10^4
% simoptions.parallel= 1; % 在并行CPU上使用稀疏矩阵求解，然后将答案传输到GPU
% vfoptions.parallel = 1;  % 1=CPU，2=GPU（默认）
% simoptions.usegpu = 0;
AgeWeightsParamNames={'mewj'}; % 许多有限视野模型对不同"年龄"应用不同权重；例如，由于生存率或人口增长率。

Params.fractionretired=sum(Params.mewj.*Params.bvec); % 注意：bvec实际上只是退休的指标

%% 设置一般均衡条件（关于资产/利率，假设具有Cobb-Douglas生产函数的代表性企业）

% 要评估的函数
FnsToEvaluate.K = @(aprime,a,z) a; % 总资产（即本期状态）
FnsToEvaluate.L = @(aprime,a,z,ybarj) exp(z+ybarj); % 总劳动供给（以效率单位计）
FnsToEvaluate.Beq = @(aprime,a,z,sj,r,tau) (1-sj)*aprime*(1+r*(1-tau)); % 总意外遗产
% 注意总劳动供给实际上是完全外生的，所以我可以预先计算它，但我感觉懒得做。

% 一般均衡方程
% 回顾GEPriceParamNames={'r','b','T'};
GeneralEqmEqns.capitalmarket = @(r,L,K,A,alpha,delta) r-(A*alpha*(K^(alpha-1))*(L^(1-alpha))-delta); % 资产回报率与资本边际产品相关
GeneralEqmEqns.SSbalance = @(b,K,L,theta,fractionretired,alpha,A) b*fractionretired-theta*(A*(1-alpha)*(K^(alpha))*(L^(-alpha)))*L; % 退休福利等于工资税收入：b*fractionretired-theta*w*L
GeneralEqmEqns.Bequests = @(T,Beq,n) T-Beq/(1+n); % 一次性转移等于意外遗产

%% 我们已经完成了模型的设置，现在我们需要为Huggett(1996)所做的所有变体求解。
% Hugget求解了多个经济体
% Params.sigma=1.5, 3
sigmavec=[1.5,3];
% Params.alowerbar=0, -w
bc_equalsminuswvec=[0,1];
% Params.sigmasqepsilon=0, 0.045
sigmasqepsilonvec=[0, 0.045];
% 确定生命周期：beta=0.994, Params.sj=1, beta=1.011, Params.sj=向量

sigma_c=1;
alowerbar_c=1;
sigmasqepsilon_c=2;
uncertainlifetime_c=2;

for sigma_c=1:2
    Params.sigma=sigmavec(sigma_c);
    for alowerbar_c=1:2
        Params.bc_equalsminusw=bc_equalsminuswvec(alowerbar_c);
        for sigmasqepsilon_c=1:2
            Params.sigmasqepsilon=sigmasqepsilonvec(sigmasqepsilon_c);
            for uncertainlifetime_c=1:2
                fprintf('Current economy is %i %i %i %i \n', sigma_c, alowerbar_c, sigmasqepsilon_c, uncertainlifetime_c)
                if uncertainlifetime_c==1 % 不确定生命周期
                    Params.beta=1.011; 
                    Params.sj=1-Params.dj;
                elseif uncertainlifetime_c==2 % 确定生命周期
                    Params.beta=0.994;
                    Params.sj=ones(size(Params.dj));
                end
                if sigma_c==1 && alowerbar_c==1 && sigmasqepsilon_c==2 && uncertainlifetime_c==2
                    heteroagentoptions.fminalgo=6; % 需要使用约束优化
                    % 强制转移必须>=0
                    heteroagentoptions.lb=[-Inf,-Inf,0]; % 回顾GEPriceParamNames={'r','b','T'};
                    heteroagentoptions.ub=[Inf,Inf,Inf];
                else
                    heteroagentoptions=struct(); % 不需要设置任何heteroagentoptions
                end
                OutputResults=Huggett1996_Fn(Params, n_a,n_z,N_j, a_grid, ReturnFn, DiscountFactorParamNames, AgeWeightsParamNames, FnsToEvaluate, GEPriceParamNames, GeneralEqmEqns,heteroagentoptions, simoptions,vfoptions,CheckUniquenessOfGE);
                FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults=OutputResults;
                counter=[sigma_c, alowerbar_c, sigmasqepsilon_c, uncertainlifetime_c];
                save ./SavedOutput/Huggett1996_Counter.mat counter
            end
        end
    end
end

save ./SavedOutput/Huggett1996_FullResults.mat FullResults
% load ./SavedOutput/Huggett1996_FullResults.mat FullResults

%% 绘制Huggett(1996)的图表

% 图1 - 收入曲线（相对于总体平均值的比率）
figure(1);
plot(19+(1:1:Params.J),exp(Params.ybarj)./sum(exp(Params.ybarj).*Params.mewj))
title({'Earnings Profile (ratio to overall mean)'})
xlabel('Age')
ylabel('Earnings')
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure1.png')

AgeConditionalStats=FullResults(1,2,2,1).OutputResults.AgeConditionalStats;
% 图2 - 不确定生命周期的财富曲线
figure(2);
plot(19+(1:1:71),AgeConditionalStats.K.Mean(1:end-8),19+(1:1:71),AgeConditionalStats.K.QuantileCutoffs(2,1:end-8),19+(1:1:71),AgeConditionalStats.K.QuantileCutoffs(5,1:end-8),19+(1:1:71),AgeConditionalStats.K.QuantileCutoffs(10,1:end-8))
legend('Mean','10th','25th','Median')
title({'Wealth Profiles: Uncertain Lifetimes'})
xlabel('Age')
ylabel('Wealth')
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure2.png')

% Huggett(1996)的图3基于美国数据。

% 图4 - 确定生命周期：年龄组内的基尼系数
figure(3);
AgeConditionalStats1=FullResults(1,1,2,2).OutputResults.AgeConditionalStats;
AgeConditionalStats2=FullResults(1,2,2,2).OutputResults.AgeConditionalStats;
plot(19+(11:1:56),AgeConditionalStats1.K.Gini(11:(end-8-15)),19+(11:1:56),AgeConditionalStats2.K.Gini(11:(end-8-15)))
legend('a=0, sigma=1.5', 'a=-w, sigma=1.5')
title({'Certain Lifetimes: Gini coefficients within age groups'})
xlabel('Age')
ylabel('Wealth Gini')
ylim([0,1.4])
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure4.png')

% 图5 - 不确定生命周期：年龄组内的基尼系数
figure(4);
AgeConditionalStats1=FullResults(1,1,2,1).OutputResults.AgeConditionalStats;
AgeConditionalStats2=FullResults(1,2,2,1).OutputResults.AgeConditionalStats;
plot(19+(11:1:56),AgeConditionalStats1.K.Gini(11:(end-8-15)),19+(11:1:61),AgeConditionalStats2.K.Gini(11:(end-8-10)))
legend('a=0, sigma=1.5', 'a=-w, sigma=1.5')
title({'Uncertain Lifetimes: Gini coefficients within age groups'})
xlabel('Age')
ylabel('Wealth Gini')
ylim([0,1.4])
saveas(gcf,'./SavedOutput/Graphs/Huggett1996_Figure5.png')


%% Huggett(1996)的表格

Table3=cell(8,9);
sigma_c=1;
alowerbarstr={'0','-w'};
% 第1行 - 确定生命周期，无收入冲击，借贷约束=0
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(1,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% 第2行 - 确定生命周期，无收入冲击，借贷约束=-w
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(2,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% 第3行 - 确定生命周期，有收入冲击，借贷约束=0
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(3,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% 第4行 - 确定生命周期，有收入冲击，借贷约束=-w
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(4,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% 第5行 - 不确定生命周期，无收入冲击，借贷约束=0
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(5,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% 第6行 - 不确定生命周期，无收入冲击，借贷约束=-w
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(6,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% 第7行 - 不确定生命周期，有收入冲击，借贷约束=0
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(7,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 8
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table3(8,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};

%Table 3
FID = fopen('./SavedOutput/LatexInputs/Huggett1996_Table3.tex', 'w');
fprintf(FID, '\\center{Wealth Distribution (risk aversion coefficient $\\sigma=$%8.1f) \\\\ \n ', sigmavec(sigma_c));
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllllll} \\hline \n');
fprintf(FID, 'Credit & Earnings &       & Transfer &         & \\multicolumn{3}{l}{Percentage wealth in}   & Zero or \\\\ \n');
fprintf(FID, 'limit  & shock    &       & wealth   & Wealth  & \\multicolumn{3}{l}{the top}                & negative \\\\ \\cline{6-8} \n');
fprintf(FID, '$\\underbar{a}$  & $\\sigma_{\\epsilon}^2$ & $K/Y$ & ration & Gini & 1\\%% & 5\\%% & 20 \\%% & wealth (\\%%) \\\\ \n');
fprintf(FID, '\\multicolumn{2}{l}{US Economy} & 3.0 & 0.78--1.32 & 0.72 & 28 & 49 & 75 & 5.8--15.0 \\\\ \n');
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Certain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{1,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{2,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{3,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{4,:});
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Uncertain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{5,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{6,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{7,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table3{8,:});
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d. \n', n_a, n_z);
fprintf(FID, '}} \\end{minipage} }');
fclose(FID);



Table4=cell(8,9);
sigma_c=2;
alowerbarstr={'0','-w'};
% Row 1
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(1,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 2
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(2,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 3
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(3,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 4
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=2;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(4,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 5
alowerbar_c=1; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(5,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 6
alowerbar_c=2; sigmasqepsilon_c=1; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(6,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 7
alowerbar_c=1; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(7,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};
% Row 8
alowerbar_c=2; sigmasqepsilon_c=2; uncertainlifetime_c=1;
OutputResults=FullResults(sigma_c,alowerbar_c,sigmasqepsilon_c,uncertainlifetime_c).OutputResults;
Table4(8,:)={alowerbarstr{alowerbar_c},sigmasqepsilonvec(sigmasqepsilon_c),OutputResults.KdivY, OutputResults.TransferWealthRatio, OutputResults.WealthGini, OutputResults.TopWealthShares(3), OutputResults.TopWealthShares(2), OutputResults.TopWealthShares(1), OutputResults.FractionWithZeroOrNegAssets};

%Table 4
FID = fopen('./SavedOutput/LatexInputs/Huggett1996_Table4.tex', 'w');
fprintf(FID, '\\center{Wealth Distribution (risk aversion coefficient $\\sigma=$%8.1f) \\\\ \n ', sigmavec(sigma_c));
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllllllll} \\hline \n');
fprintf(FID, 'Credit & Earnings &       & Transfer &         & \\multicolumn{3}{l}{Percentage wealth in}   & Zero or \\\\ \n');
fprintf(FID, 'limit  & shock    &       & wealth   & Wealth  & \\multicolumn{3}{l}{the top}                & negative \\\\ \\cline{6-8} \n');
fprintf(FID, '$\\underbar{a}$  & $\\sigma_{\\epsilon}^2$ & $K/Y$ & ration & Gini & 1\\%% & 5\\%% & 20 \\%% & wealth (\\%%) \\\\ \n');
fprintf(FID, '\\multicolumn{2}{l}{US Economy} & 3.0 & 0.78--1.32 & 0.72 & 28 & 49 & 75 & 5.8--15.0 \\\\ \n');
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Certain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{1,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{2,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{3,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{4,:});
fprintf(FID, '\\multicolumn{9}{l}{\\textit{Uncertain Lifetimes}} \\\\ \n');
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{5,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{6,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{7,:});
fprintf(FID, '%s & %8.3f & %8.1f & %8.2f & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', Table4{8,:});
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Based on baseline grid size for assets of %d, and shocks of %d. \n', n_a, n_z);
fprintf(FID, '}} \\end{minipage} }');
fclose(FID);





























