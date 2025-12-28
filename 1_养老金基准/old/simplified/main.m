clear
warning off
% 主程序

% 设置参数
tb = 20;     % initial starting age
td = 100;    % death age
tn = td-tb+1;% total period
na = 51;
ncash = 51;
n = 5;

maxcash = 50.0;
mincash = 0.25;

% 基础收入系数
aa = -2.170042+2.700381;
b1 = 0.16818;
b2 = -0.0323371/10;
b3 = 0.0019704/100;

% 其他参数设置
delta = 0.97;
r = 1.015;
mu = 0.04;
sigr = 0.2;
riskaversion = 3.84;

% 生成生存概率
survprob = [ ...
    0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, ...
    0.99963, 0.99961, 0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, ...
    0.99934, 0.99928, 0.99922, 0.99916, 0.99908, 0.99900, 0.99891, 0.99881, ...
    0.99870, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, 0.99776, 0.99756, ...
    0.99735, 0.99713, 0.99690, 0.99666, 0.99640, 0.99613, 0.99583, 0.99551, ...
    0.99515, 0.99476, 0.99432, 0.99383, 0.99330, 0.99270, 0.99205, 0.99133, ...
    0.99053, 0.98961, 0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, ...
    0.97391, 0.96943, 0.96429, 0.95854, 0.95221, 0.94537, 0.93805, 0.93027, ...
    0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, 0.85846, 0.84452, ...
    0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.73240, 0.70893, 0.68424 ...
    ]';

% 生成网格点
[grid, weig] = tauchen(n, 0, 0, 1, 2);
weig = diag(weig);

% 计算风险资产收益
gret = zeros(n,1);
for i1 = 1:n
    gret(i1) = r + mu + grid(i1)*sigr;
end

% 生成cash-on-hand网格
l_maxcash = log(maxcash);
l_mincash = log(mincash);
stepcash = (l_maxcash-l_mincash)/(ncash-1);

lgcash = zeros(ncash,1);
gcash = zeros(ncash,1);
for i1 = 1:ncash
    lgcash(i1) = l_mincash+(i1-1)*stepcash;
    gcash(i1) = exp(lgcash(i1));
end

% 初始化策略函数和值函数
C = zeros(ncash,tn);
A = ones(ncash,tn);
V = zeros(ncash,tn);

% 最后一期的值函数
C(:,end) = 1;
A(:,end) = 0;
for i1 = 1:ncash
    V(i1,end) = (C(i1,end)*gcash(i1))^(1-riskaversion)/(1-riskaversion);
end

% 计算永久收入冲击
yh = zeros(n, 1);
for i = 1:n
    yh(i) = exp(grid(i));
end

% 向后递归求解
for i1 = 2:tn-1
    t = tn - i1 + 1;
    fprintf('Solving period %d\n', t);

    % 计算当期收入
    y_base = f_y(t);

    % 创建收入和收益率的组合
    [Y, R] = meshgrid(yh, gret);
    nweig1 = kron(weig, ones(size(Y)));

    for i3 = 1:ncash
        x0 = [0.5; 0.2];
        lb = [0; 0];
        ub = [1; 1];

        % 使用fmincon优化
        options = optimoptions('fmincon','Display','off');
        [x,fval] = fmincon(@(x)cocco_fun_valuefunc(x,gcash(i3),V(:,t),gret,r,y_base+yh,delta,gcash,survprob(t-1),...
            nweig1,riskaversion,n),x0,[],[],[],[],lb,ub,[],options);

        C(i3,t-1) = x(1);
        if x(1) == 1
            A(i3,t-1) = NaN;
        else
            A(i3,t-1) = x(2);
        end
        V(i3,t-1) = -fval;
    end
end

% 保存结果
save('cocco.mat','A','V','C');

% Tauchen方法离散化AR(1)过程
function [Z, Zprob] = tauchen(N, mu, rho, sigma, m)
Z = zeros(N, 1);
Zprob = zeros(N, N);
a = (1 - rho) * mu;

Z(end) = m * sqrt(sigma^2 / (1 - rho^2));
Z(1) = -Z(end);
zstep = (Z(end) - Z(1)) / (N - 1);

for i = 2:N
    Z(i) = Z(1) + zstep * (i-1);
end

Z = Z + a / (1 - rho);

for j = 1:N
    for k = 1:N
        if k == 1
            Zprob(j,k) = normcdf((Z(1) - a - rho * Z(j) + zstep/2) / sigma);
        elseif k == N
            Zprob(j,k) = 1 - normcdf((Z(end) - a - rho * Z(j) - zstep/2) / sigma);
        else
            up = normcdf((Z(k) - a - rho * Z(j) + zstep/2) / sigma);
            down = normcdf((Z(k) - a - rho * Z(j) - zstep/2) / sigma);
            Zprob(j,k) = up - down;
        end
    end
end
end


% 计算基础收入函数
function y = f_y(t)
% 基础收入系数
aa = -2.170042 + 2.700381;
b1 = 0.16818;
b2 = -0.0323371/10;
b3 = 0.0019704/100;

% 计算基础收入
y = aa + b1*t + b2*t^2 + b3*t^3;
end

% 值函数计算
function v = cocco_fun_valuefunc(x, cash, nextV, gret, rf, income, delta, gcash, survprob, weig, riskaversion, n)
    auxVV = 0;
    u = (x(1)*cash)^(1-riskaversion)/(1-riskaversion);
    
    % 循环计算版本
    for i6 = 1:n
        for i7 = 1:n
            % 下期的cash-on-hand
            sav = (cash*(1 - x(1)));
            cash_1 = (rf * (1 - x(2)) + gret(i6) * x(2)) * sav + income(i7);
            cash_1 = min(max(cash_1, gcash(1)), gcash(end)); % 限制在gcash的范围内
            % 插值计算
            int_V = interp1(gcash(:), nextV, cash_1,'spline');
            auxVV = auxVV + weig(i6, i7) * survprob * int_V;
        end
    end

    v = -(u + delta * auxVV);
end