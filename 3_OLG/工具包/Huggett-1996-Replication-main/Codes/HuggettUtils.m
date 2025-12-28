classdef HuggettUtils
    methods (Static)
        % HHPrices_Huggett函数
        function [Y, R, w, b] = HHPrices_Huggett(K, L, cS)
            % 生产函数相关计算
            Y   = cS.A * (K^cS.alpha) * (L^(1-cS.alpha));
            MPK = cS.alpha * cS.A * (K^(cS.alpha-1)) * (L^(1-cS.alpha));
            MPL = (1-cS.alpha) * cS.A * (K^cS.alpha) * (L^(-cS.alpha));

            % 税收
            tau = 0.195/(1-cS.ddk * K/Y);

            % 家户面临的价格
            R   = 1 + (MPK - cS.ddk)*(1 - tau);
            w   = MPL*(1 - tau - cS.theta);

            % 社会保障福利
            b   = cS.theta * w * L/cS.retireMass;
        end
        
        % ParameterValues_Fixed函数
        function cS = ParameterValues_Fixed()
            %% 人口统计学参数
            % 物理年龄
            cS.age1      = 20;     % 模型经济中的每个人从20岁开始
            cS.ageLast   = 98;     % 最大年龄
            cS.ageRetire = 65;     % 家户在65岁退休
            cS.popGrowth = 0.012;  % 人口增长率

            % 模型年龄
            cS.aD        = cS.ageLast - cS.age1 + 1;
            cS.aR        = cS.ageRetire - cS.age1 + 1;

            % 生存概率
            % 基于Jordan, C., Life contingencies, 2nd ed. (Society of Actuaries)
            % 条件死亡概率：物理年龄20到98岁
            cS.d         = [0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
                            0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
                            0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
                            0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
                            0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
                            0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
                            0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
                            0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018]; 
            % 条件生存概率
            cS.s         = 1 - cS.d; 

            % 按年龄划分的家户数量
            cS.ageMassV   = ones(1, cS.aD);
            for i = 2 : length(cS.ageMassV)
                cS.ageMassV(i) = cS.s(i-1) * cS.ageMassV(i-1) / (1 + cS.popGrowth);
            end
            cS.ageMassV   = cS.ageMassV./sum(cS.ageMassV);

            % 退休家户数量
            cS.retireMass = sum(cS.ageMassV(cS.aR + 1 : end));

            % 每个模型年龄的物理年龄
            cS.physAgeV   = (cS.age1 : cS.ageLast)';

            %% 家户参数
            cS.sigma      = 1.5;   % 效用函数曲率
            cS.beta       = 1.011; % 折现因子
            cS.cFloor     = 0.05;  % 消费下限（数值原因引入）
            cS.nSim       = 5e4;   % 模拟个体数量

            %% 技术参数
            cS.A          = 0.895944;
            cS.alpha      = 0.36;
            cS.ddk        = 0.06;

            %% 社会保障参数
            cS.theta      = 0.1;

            %% 劳动禀赋参数
            cS.leSigma1      = 0.38 ^ 0.5;
            cS.leShockStd    = 0.045 .^ 0.5;
            cS.lePersistence = 0.96;     
            cS.leWidth       = 4;     % 网格宽度（标准差数量）
            cS.nw            = 18;    % 劳动禀赋网格大小

            %% 网格参数
            % 资本网格
            cS.tgKY          = 3;          % 目标资本/产出比为3
            cS.tgWage        = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk            = 50;
            cS.kMin          = 0;          % 由于借贷约束：k>=0
            cS.kMax 	     = 100 * cS.tgWage;
            kGridV           = linspace(cS.kMin, cS.kMax, cS.nk); 
            cS.kGridV        = kGridV(:);  % 转为列向量
        end
        
        % tauchen函数
        function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std)
            % 使用网格宽度
            a_bar = n_std * sqrt(pSigma^2 / (1 - pRho^2));

            % 网格
            y     = linspace(-a_bar, a_bar, N);

            % 点之间的距离
            d     = y(2) - y(1);

            % 获取转移概率
            trProbM  = zeros(N, N);

            for iRow = 1 : N
               % 先处理端点
               trProbM(iRow, 1) = normcdf((y(1) - pRho*y(iRow) + d/2) / pSigma);
               trProbM(iRow, N) = 1 - normcdf((y(N) - pRho*y(iRow) - d/2) / pSigma);

               % 填充中间列
               for iCol = 2 : N-1
                  trProbM(iRow, iCol) = (normcdf((y(iCol) - pRho*y(iRow) + d/2) / pSigma) - ...
                                         normcdf((y(iCol) - pRho*y(iRow) - d/2) / pSigma));
               end
            end

            % 围绕均值（ybar / (1 - rho)）中心化过程
            y = y + pMu / (1 - pRho); 
        end
        
        % norm_grid函数
        function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig)
            n = length(xV);
            
            % 验证输入
            if any(xV(2:n) < xV(1:n-1))
                error('xV not increasing');
            end
            
            if xMin > xV(1) || xMax < xV(n)
                error('Invalid xMin or xMax');
            end
            
            if mu < xMin || mu > xMax
                error('Invalid mu');
            end

            % 构建区间边界
            xMidV = 0.5 .* (xV(1:n-1) + xV(2:n));
            lbV   = [xMin; xMidV(:)];
            ubV   = [xMidV(:); xMax];

            % 找出每个区间的质量
            cdfV  = normcdf([xMin; ubV], mu, sig);
            massV = cdfV(2:(n+1)) - cdfV(1:n);
            massV = massV(:) ./ sum(massV);
            
            % 输出检查
            if any(ubV < lbV)
                error('ubV < lbV');
            end
        end
        
        % LaborSupply_Huggett函数
        function [HHlaborM, L] = LaborSupply_Huggett(eIdxM, cS, paramS)
            % 个体劳动供给：效率 * 劳动禀赋冲击
            HHlaborM = zeros(cS.nSim, cS.aD);
            for a = 1 : cS.aD
               HHlaborM(:, a) = paramS.ageEffV(a) .* paramS.leGridV(eIdxM(:,a));
            end

            % 总劳动供给
            L = mean(HHlaborM,1)* cS.ageMassV';
        end
        
        % EarningProcess_olgm函数
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            [logGridV, trProbM] = HuggettUtils.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);

            % 新代理从近似对数正态分布中抽取
            prob1V = HuggettUtils.norm_grid(logGridV, logGridV(1)-2, logGridV(end)+2, 0, cS.leSigma1);
            prob1V = prob1V(:);

            % 改进缩放
            logGridV = logGridV(:) - logGridV(1) - 1;
            
            % 输出验证
            if any(abs(sum(trProbM, 2) - 1) > 1e-6)
                error('probs do not sum to 1');
            end
            
            if abs(sum(prob1V) - 1) > 1e-6
                error('prob1V does not sum to 1');
            end
        end
        
        % MarkovChainSimulation函数
        function eIdxM = MarkovChainSimulation(nSim, T, prob0V, trProbM, rvInM)
            ns = length(prob0V);
            
            % 输入验证
            if max(abs(sum(trProbM) - 1)) > 1e-5
                error('Probabilities do not sum to one');
            end
            
            if abs(sum(prob0V) - 1) > 1e-5
                error('Initial probs do not sum to 1');
            end
                
            % 准备工作
            cumTrProbM = cumsum(trProbM);
            cumTrProbM(ns, :) = 1;
            cumTrProbM = cumTrProbM';
            
            % 按日期迭代
            eIdxM = zeros([nSim, T]);
            
            % 抽取t=1
            eIdxM(:, 1) = 1 + sum((rvInM(:,1) * ones(1, ns)) > (ones(nSim,1) * cumsum(prob0V(:)')), 2);
            
            % 对于t = 2, ..., T
            for t = 1 : (T-1)
               eIdxM(:, t+1) = 1 + sum((rvInM(:,t+1) * ones(1, ns)) > cumTrProbM(eIdxM(:,t), :), 2);
            end
        end
        
        % LaborEndowSimulation_olgm函数
        function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
            % 设置随机数生成器以便重复
            rng(433);
            
            % 按[ind, age]的禀赋状态
            eIdxM = HuggettUtils.MarkovChainSimulation(cS.nSim, cS.aD, paramS.leProb1V, ...
                                      paramS.leTrProbM', rand([cS.nSim, cS.aD]));
        end
        
        % HHSimulation_olgm函数
        function [kHistM, cHistM] = HHSimulation_olgm(kPolM, cPolM, eIdxM, cS)
            % 输入验证
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               'size', [cS.nk,cS.nw,cS.aD]});
            
            % 模拟资本和消费历史，按年龄
            nSim   = size(eIdxM, 1);
            kHistM = zeros(nSim, cS.aD);
            cHistM = zeros(nSim, cS.aD);
            
            for a = 1 : cS.aD
                for ie = 1 : cS.nw
                  % 找到此年龄劳动禀赋为ie的家户
                  idxV = find(eIdxM(:,a) == ie);
            
                  if ~isempty(idxV)
                     if a < cS.aD
                        % 通过插值找到每个个体的下一期资本
                        kHistM(idxV, a+1) = interp1(cS.kGridV(:), kPolM(:,ie,a), ...
                                                    kHistM(idxV, a), 'linear');
                     end
                     
                     cHistM(idxV, a) = interp1(cS.kGridV(:), cPolM(:,ie,a), ...
                                               kHistM(idxV, a), 'linear');
                  end
            
                end
            end
            
            % 输出验证
            validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               '>', cS.kGridV(1) - 1e-6});
        end
        
        % CES_utility函数
        function [muM, utilM] = CES_utility(cM, sig)
            % 输入验证
            if any(cM(:) < 1e-8)
                error('Cannot compute utility for very small consumption');
            end
            
            if length(sig) ~= 1
                error('sig must be scalar');
            end
            
            if sig <= 0
                error('sig must be > 0');
            end
            
            % 计算效用
            if sig == 1                            % 对数效用
               utilM = log(cM);                    % 效用
               muM   = 1 ./ cM;                    % 边际效用
            else
               utilM = cM .^ (1-sig) ./ (1-sig);   % CES效用
               muM = cM .^ (-sig);                 % 边际效用
            end
        end
        
        % HHIncome_Huggett函数
        function incomeM = HHIncome_Huggett(kV, R, w, T, b, a, paramS, cS)
            % 非资本收入（按冲击）
            if a <= cS.aR
               nonCapIncomeV = paramS.ageEffV(a) .* w .* paramS.leGridV + T;
            else
               nonCapIncomeV = b .* ones([cS.nw, 1]) + T;
            end
            
            % 给定年龄的总收入
            incomeM = R * kV(:) * ones([1, cS.nw]) + ones([length(kV),1]) * nonCapIncomeV(:)';
        end
        
        % HHSolution_VFI_Huggett函数
        function [cPolM, kPolM, valueM] = HHSolution_VFI_Huggett(R, w, T, bV, paramS, cS)
            % 初始化策略函数
            cPolM  = zeros(cS.nk, cS.nw, cS.aD);
            kPolM  = zeros(cS.nk, cS.nw, cS.aD);
            valueM = zeros(cS.nk, cS.nw, cS.aD);
            
            % 向后归纳
            for a = cS.aD : -1 : 1
               % 下一期值函数
               if a < cS.aD
                  vPrimeM = valueM(:,:,a+1);
               else
               % 没有下一期，因为a=aD是最后一期
                  vPrimeM = [];
               end
            
               [cPolM(:,:,a), kPolM(:,:,a), valueM(:,:,a)] = ...
                  HuggettUtils.HHSolutionByAge_VFI_Huggett(a, vPrimeM, R, w, T, bV(a), paramS, cS);
              
            end
            
            % 输出验证
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'positive', 'size', [cS.nk, cS.nw, cS.aD]});
            
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'size', [cS.nk, cS.nw, cS.aD]});
            
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'size', [cS.nk, cS.nw, cS.aD], ...
                               '>=', cS.kMin - 1e-6, '<=', cS.kGridV(end) + 1e-6});
        end
        
        % HHSolutionByAge_VFI_Huggett函数
        function [cPolM, kPolM, valueM] = HHSolutionByAge_VFI_Huggett(a, vPrime_keM, R, w, T, b, paramS, cS)
            % 输入检查
            if a < cS.aD
                if ~isequal(size(vPrime_keM), [cS.nk, cS.nw])
                    error('Invalid size of cPrimeM');
                end
            end
            
            % 年龄a时的收入y
            yM = HuggettUtils.HHIncome_Huggett(cS.kGridV, R, w, T, b, a, paramS, cS);
            
            % 优化选项
            fminbndOptS = optimset('fminbnd');
            fminbndOptS.TolX = 1e-5;
            
            if a == cS.aD
               % 吃掉所有收入，不储蓄，因为这是最后一期
               cPolM = yM;
               kPolM = zeros(cS.nk, cS.nw);
               [~, valueM] = HuggettUtils.CES_utility(cPolM, cS.sigma);
               
            else
               % 为策略函数和值函数分配空间
               cPolM = zeros(cS.nk, cS.nw);
               kPolM = zeros(cS.nk, cS.nw);
               valueM = zeros(cS.nk, cS.nw);
            
               % 遍历状态[ik, ie]
               for ie = 1 : cS.nw
                   % 期望值函数，按kPrime网格点 -- EV(k')
                   ExValuePrimeV = zeros(cS.nk, 1);
                   for ik = 1 : cS.nk
                      ExValuePrimeV(ik) = paramS.leTrProbM(ie,:) * vPrime_keM(ik,:)';
                   end
            
                   % 明天EV(k')的连续近似
                   vPrimeOfK = griddedInterpolant(cS.kGridV, ExValuePrimeV, 'linear');
            
                   % 遍历资本状态
                   for ik = 1 : cS.nk
                        [cPolM(ik,ie), kPolM(ik,ie), valueM(ik,ie)] = ...
                              HuggettUtils.HHSolutionByOneState_VFI_Huggett(a, yM(ik,ie), R, vPrimeOfK, fminbndOptS, cS);
                   end
            
               end
            end
            
            % 输出验证
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'positive', 'size', [cS.nk, cS.nw]});
            
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               '>=', cS.kMin - 1e-6,  '<=', cS.kGridV(cS.nk) + 1e-6, ...
                               'size', [cS.nk, cS.nw]});
            
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'size', [cS.nk, cS.nw]});
        end
        
        % HHSolutionByOneState_VFI_Huggett函数
        function [c, kPrime, ValueFunc] = HHSolutionByOneState_VFI_Huggett(a, y, R, vPrimeOfK, fminbndOptS, cS)
            % 可行kPrime的范围
            kPrimeMax = min(cS.kGridV(cS.nk), y - cS.cFloor);
            
            if kPrimeMax <= cS.kGridV(1)
               % 没有可行的选择。家户获得最低消费且不储蓄
               kPrime = cS.kGridV(1);
               
            else
               % 找到最优kPrime
               [kPrime, ~, ~] = fminbnd(@Bellman, cS.kGridV(1), kPrimeMax, fminbndOptS);
               
            end
            
            [ValueFunc, c] = Bellman(kPrime);
            ValueFunc = -ValueFunc;
            
            % 输出验证
            validateattributes(kPrime, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar', ...
                               '>=', cS.kGridV(1) - 1e-6, '<=', kPrimeMax + 1e-6});
            
            validateattributes(c, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar', ...
                               '>=', cS.cFloor - 1e-6});
            
            validateattributes(ValueFunc, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar'});
                
            % 嵌套：目标函数
            function [Valfunc, c] = Bellman(kPrime)
                c = max(cS.cFloor, y - kPrime);
                [~, u] = HuggettUtils.CES_utility(c, cS.sigma);
                Valfunc = -(u + cS.beta * cS.s(a) * R * vPrimeOfK(kPrime));
            end
        end
    end
end 