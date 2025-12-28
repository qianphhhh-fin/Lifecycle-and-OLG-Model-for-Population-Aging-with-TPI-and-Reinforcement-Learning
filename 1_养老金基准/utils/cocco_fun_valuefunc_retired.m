function v = cocco_fun_valuefunc_retired(x,cash,nextV,gret,rf,ret_fac,gamma,beta,psi_1,psi_2,theta,gcash,survprob,weig)
% ============= 输入 ===============
% ==== 决策变量：====
% x(1)-----消费
% x(2)-----风险资产投资比例
% x(3)-----养老金购买比例
% ==== 状态变量: ====
% cash-----手中现金
% fund-----养老基金余额
% 状态变量网格：gcash,gfund
% ==== 外生参数: ====
% gret-----随机收益率
% rf-------无风险收益率
% nextV----下期的值函数grid
% zeta-----养老金提取比例
% ret_fac--固定养老金
% survprob--当期存活概率
% n--------随机状态个数
% weig-----每个随机状态的权重
% ============= 输出 ===============
% v--------函数值

auxVV = 0;
sav = cash-x(1);
delta = 0.2;
% u=(x(1))^(1-gamma)/(1-gamma);
% if x(1)<delta
%     C = -gamma/2*delta^(-gamma-1);
%     B = delta^(-gamma)-2*C*delta;
%     A = delta^(1-gamma)/(1-gamma)-B*delta-C*delta^2;
%     u = A+B*x(1)+C*x(1)^2;
% else
u = -(x(1))^(1-gamma);
% end
 % u=log(x(1));
% if u<-2
%     u =-2 ;
% end
% u=exp((x(1))^(1-gamma)/(1-gamma));

% ====== 下期的cash-on-hand ====================
% loc = [];
% 非向量化操作
cash_1 = (rf*(1-x(2))+gret.*x(2))*sav + ret_fac ; % 下期的cash-on-hand
% ====== 下期的pension fund ====================
% fund_1 = fund; % 养老金不再计息

% 插值
% permute用作调换V的维度顺序以符合interp3的要求
% Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
%  X、Y 和 Z 是网格向量，则 size(V) = [length(Y) length(X) length(Z)]。
% V 原本的顺序[cash,fund,house] -> 新的顺序[fund cash house]
int_V  = interp1(gcash,nextV,cash_1,'spline');
% int_V(loc) = -999;
% 不允许负的cash-on-hand
% auxVV = auxVV + weig'*survprob*(int_V.^(1-rho)); % 累加auxVV得到EV
auxVV = auxVV + weig'*survprob*(int_V); % 累加auxVV得到EV
% if auxVV == 0
%     v = -((u).^psi_2); % 
% else
%     v = -((u + beta.*(auxVV.^(1/theta))).^psi_2); % 
% end
if auxVV == 0
    v = -((u)); % 
else
    v = -((u + beta.*(auxVV))); % 
end

end
