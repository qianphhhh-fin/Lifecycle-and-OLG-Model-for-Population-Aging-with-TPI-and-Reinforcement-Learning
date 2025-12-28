function OutputResults=Huggett1996_Fn(Params, n_a,n_z,N_j, a_grid, ReturnFn, DiscountFactorParamNames, AgeWeightsParamNames, FnsToEvaluate, GEPriceParamNames, GeneralEqmEqns, heteroagentoptions, simoptions, vfoptions,CheckUniquenessOfGE)
% Huggett (1996) 生命周期经济中的财富分配模型复现

%% 需要创建适当的 z_grid 和 pi_z
Params.sigmasqy=Params.sigmasqepsilon./(1-Params.gamma.^2); 
% 收入的初始分布
Params.sigmasqy1=0.38;

% Huggett 使用17个在 +-4*sigmasqy1 之间均匀分布的状态，外加第18个
% 状态在 6*sigmasqy1，"状态之间的转移概率是通过
% 对当前状态值条件下的正态分布下方区域进行积分计算的。"
z_grid=[linspace(-4*sqrt(Params.sigmasqy1),4*sqrt(Params.sigmasqy1),17),6*sqrt(Params.sigmasqy1)]';
pi_z=nan(18,18);
% 以下几行实现转移矩阵，主要是从 TauchenMethod() 命令复制的代码
sigma=sqrt(Params.sigmasqepsilon); %e的标准差
for ii=1:length(z_grid)
    pi_z(ii,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2-Params.gamma*z_grid(ii),0,sigma);
    for jj=2:(length(z_grid)-1)
        pi_z(ii,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2-Params.gamma*z_grid(ii),0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2-Params.gamma*z_grid(ii),0,sigma);
    end
    pi_z(ii,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2-Params.gamma*z_grid(ii),0,sigma);
end
% z_grid=gpuArray(z_grid);
% pi_z=gpuArray(pi_z);

%% 出生时(j=1)个体的初始分布(出生个体的初始资产为零，个体受到18种冲击的概率（行的和为1）
% jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist=zeros(n_a,n_z);
[trash,zeroassets_index]=min(abs(a_grid));
% 基于 Params.sigmasqy1（年龄1时收入的方差）
% 在正态分布假设下将它们放在现有的 z_grid 上
% 以下几行实现这一点，主要是从 TauchenMethod() 命令复制的代码
sigma=sqrt(Params.sigmasqy1); %e的标准差
for ii=1:length(z_grid)
    jequaloneDist(zeroassets_index,1)=normcdf(z_grid(1)+(z_grid(2)-z_grid(1))/2,0,sigma);
    for jj=2:(length(z_grid)-1)
        jequaloneDist(zeroassets_index,jj)=normcdf(z_grid(jj)+(z_grid(jj+1)-z_grid(jj))/2,0,sigma)-normcdf(z_grid(jj)-(z_grid(jj)-z_grid(jj-1))/2,0,sigma);
    end
    jequaloneDist(zeroassets_index,end)=1-normcdf(z_grid(end)-(z_grid(end)-z_grid(end-1))/2,0,sigma);
end

%% 求解一般均衡
% 使用工具包找到均衡价格指数。有两种方法
% 可以做到这一点。在下面，我使用"搜索"方法计算
% （初始）一般均衡。使用搜索后，我会在 p_grid 周围
% 进行一系列探索，以"确定"这是唯一的
% 均衡，或者至少是附近唯一的均衡。

% 不使用 p_grid，只是搜索。使用 n_p=0。（设置实际算法
% 为"搜索"可以通过 heteroagentoptions.fminalgo 完成）
heteroagentoptions.verbose=1;
[p_eqm,~, GeneralEqmEqnsValues]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, 0, pi_z, 0, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm.r;
Params.b=p_eqm.b;
Params.T=p_eqm.T;
fprintf('搜索找到均衡 r=%8.2f, b=%8.2f, T=%8.2f \n')
save ./SavedOutput/Huggett1996.mat Params
% load ./SavedOutput/Huggett1996.mat Params

if CheckUniquenessOfGE==1
    % 使用 p_grid。如果你想查找多重均衡的可能性，这会很有帮助。
    % 我们将其用作之前找到的均衡的非正式二次检查
    % GEPriceParamNames={'r','b','T'}; % 已在上面声明
    r_grid=linspace(0.5,2,21)'*Params.r; %可以
    b_grid=linspace(0.5,2,21)'*Params.b; %可以
    T_grid=linspace(0.5,2,21)'*Params.T; %好的
    p_grid=[r_grid,b_grid, T_grid];
    %
    disp('计算对应于稳态均衡的价格向量')
    n_p=[length(r_grid),length(b_grid),length(T_grid)];
    heteroagentoptions.pgrid=p_grid;
    heteroagentoptions.verbose=1;
    [p_eqm,p_eqm_index1, GeneralEqmEqnsValues1]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, n_p, pi_z, 0, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
    Params.r=p_eqm.r;
    Params.b=p_eqm.b;
    Params.T=p_eqm.T;
    fprintf('粗网格搜索找到均衡 r=%8.2f, b=%8.2f, T=%8.2f \n')
    
    % 再次网格搜索，只是为了更加确定和准确
    r_grid=linspace(0.95,1.05,21)'*Params.r; %可以
    b_grid=linspace(0.95,1.05,21)'*Params.b; %可以
    T_grid=linspace(0.95,1.05,21)'*Params.T; %好的
    p_grid=[r_grid,b_grid, T_grid];
    [p_eqm,p_eqm_index2, GeneralEqmEqnsValues2]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightsParamNames,0, n_a, n_z, N_j, n_p, pi_z, 0, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
    Params.r=p_eqm.r;
    Params.b=p_eqm.b;
    Params.T=p_eqm.T;
    fprintf('细网格搜索找到均衡 r=%8.2f, b=%8.2f, T=%8.2f \n')
end

%% 计算均衡的一些特性
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,0,n_a,n_z,N_j,pi_z,Params,simoptions);

% 从一些基础开始（这些之前在返回函数内部完成）：
% 重新排列 r=MPK-delta 得到以下方程（MPK是资本边际产出）
KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% 从K/Y，代入生产函数得到
KdivY=(KdivL^(1-Params.alpha))/Params.A;
% 我们知道 w=MPL（MPL是劳动边际产出）
Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % 工资率（每单位有效劳动）
% Huggett (1996) 将tau校准为以下值（见第478页的解释）
Params.tau=0.195/(1-Params.delta*KdivY);

% 总财富转移
% 从一些基础开始（这些之前在返回函数内部完成）：
% 重新排列 r=MPK-delta 得到以下方程（MPK是资本边际产出）
KdivL=((Params.r+Params.delta)/(Params.alpha*Params.A))^(1/(Params.alpha-1));
% 从K/Y，代入生产函数得到
KdivY=(KdivL^(1-Params.alpha))/Params.A;
% 我们知道 w=MPL（MPL是劳动边际产出）
Params.w=Params.A*(1-Params.alpha)*(KdivL^Params.alpha); % 工资率（每单位有效劳动）
% Huggett (1996) 将tau校准为以下值（见第478页的解释）
Params.tau=0.195/(1-Params.delta*KdivY);

% 总财富转移
AggregateWealthTransers=zeros(1,N_j);
for jj=1:Params.J
    for ii=1:jj
        AggregateWealthTransers(jj)=AggregateWealthTransers(jj)+Params.T*(1+Params.r*(1-Params.tau))^(ii-1);
    end
end
AggregateWealthTransers=sum(Params.mewj.*AggregateWealthTransers);
% 总财富
FnsToEvaluate2.TotalWealth = @(aprime_val,a_val,z_val) a_val;
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate2, Params, [], 0, n_a, n_z,N_j, 0, a_grid, z_grid);
% 转移财富比率 - 这对应于Kotlikoff和Summers(1981)的概念，衡量代际转移在总财富中的比例
TransferWealthRatio=AggregateWealthTransers/AggVars.TotalWealth.Mean;


% 计算零或负财富人口比例 - 这是Huggett论文中讨论的重要指标之一
FnsToEvaluate3.ZeroOrNegAssets = @(aprime_val,a_val,z_val) (a_val<=0); % 零或负资产的指标
FractionWithZeroOrNegAssets=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate3, Params, [], 0, n_a, n_z,N_j, 0, a_grid, z_grid);

% 计算财富洛伦兹曲线（从而得到所有百分位份额）以及
% 收入的洛伦兹曲线（在第480页底部，第481页顶部的文本中，
% 有一系列关于模型收入的描述，以工作年龄为条件）
AllLorenzCurves=EvalFnOnAgentDist_LorenzCurve_FHorz_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], 0, n_a, n_z,N_j, 0, a_grid, z_grid);
TopWealthShares=100*(1-AllLorenzCurves.K([80,95,99],1)); % 需要前20%、5%和1%的份额用于Huggett(1996)的表格
% 计算财富基尼系数 - 基尼系数是洛伦兹曲线与45度线之间的面积占45度线以下面积的百分比
WealthGini=Gini_from_LorenzCurve(AllLorenzCurves.K(:,1));

AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,0,n_a,n_z,N_j,0,a_grid,z_grid);

%% 整理我想作为输出返回的一系列内容

OutputResults.Params=Params;

OutputResults.AgeConditionalStats=AgeConditionalStats;

OutputResults.WealthGini=WealthGini;
OutputResults.TopWealthShares=TopWealthShares;
OutputResults.FractionWithZeroOrNegAssets=100*FractionWithZeroOrNegAssets.ZeroOrNegAssets.Mean;
OutputResults.TransferWealthRatio=TransferWealthRatio;
OutputResults.KdivY=KdivY;

OutputResults.z_grid=gather(z_grid);
OutputResults.pi_z=gather(pi_z);

OutputResults.GeneralEqmEqnsValues=gather(GeneralEqmEqnsValues); % 用于检查均衡的准确性

if CheckUniquenessOfGE==1
    % 一些用于查看均衡是否唯一的内容
    OutputResults.p_eqm_index1=p_eqm_index1;
    OutputResults.GeneralEqmEqnsValues1=GeneralEqmEqnsValues1;
    OutputResults.p_eqm_index2=p_eqm_index2;
    OutputResults.GeneralEqmEqnsValues2=GeneralEqmEqnsValues2;
end

OutputResults=gather(OutputResults); % 确保所有内容都在cpu上，而不是gpu上



end























