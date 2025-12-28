function v = cocco_fun_valuefunc_work(x,cash,gyp,nextV,yh,gret,rf,gamma,beta,psi_1,psi_2,theta,gcash,survprob,n,weig)
% ============= 输入 ===============
% ==== 决策变量：====
% x(1)-----消费
% x(2)-----风险资产投资比例
% x(3)-----养老金购买比例
% ==== 状态变量: ====
% cash-----手中现金
% fund-----养老基金余额
% gyp------
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

 
 % u=log(x(1));% if u<-2
%     u =-2 ;
% end
% u=exp((x(1))^(1-gamma)/(1-gamma));

% u=x(1).^(1-psi_1)./(1-psi_1);
for i6=1:n % z_t, 来自P_t=P_t-1+z_t
    for i8=1:n % epsilon_t, 即风险收益随机项
        for i7=1:n % u_t, 来自logY_t = f_t + P_t + u_t
            % ====== 下期的cash-on-hand ====================
            sav = (cash-x(1))/gyp(i6,i8);
            cash_1 = (rf*(1-x(2))+gret(i8)*x(2))*sav + yh(i7,i8); % 下期的cash-on-hand
            cash_1 = max(min(cash_1,gcash(end,1)),gcash(1,1));
            % 插值
            int_V  = interp1(gcash,nextV,cash_1,'spline')';
            % auxVV = auxVV + weig(i6,i7,i8)*survprob*((int_V*gyp(i6,i8))^(1-rho)); % 累加auxVV得到EV
            auxVV = auxVV + weig(i6,i7,i8)*survprob*((int_V*gyp(i6,i8))); % 累加auxVV得到EV
        end
    end
end
% v = -(((1-beta)*u + beta.*(auxVV.^(1/theta))).^psi_2); %  
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
