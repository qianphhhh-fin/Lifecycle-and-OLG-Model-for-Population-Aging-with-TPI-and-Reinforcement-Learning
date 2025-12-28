% --- START OF FILE main_olg_v1_utils.m ---

classdef main_olg_v1_utils
    % 扩展OLG模型工具类 (v5 - Full Model with Debug Flag)
    % 针对 16 期模型，加入大量调试信息, 由 debugMode 控制

    methods (Static)
        function cS = initParameters()
            % --- Debug Flag ---
            cS.debugMode = 0; % Set to 1 to enable debug prints and saves, 0 to disable
            % ------------------

            if cS.debugMode == 1, fprintf('[DEBUG] Initializing Parameters...\n'); end
            % 人口参数
            cS.aD = 16; cS.retirementAgeIdx = 9; cS.workingAgeMaxIdx = 8;
            cS.age1 = 22; cS.ageLast = 101;
            cS.beta = [0.995, 0.99, 0.985, 0.98, 0.975, 0.97, 0.965, 0.96, ...
                       0.95, 0.94, 0.92, 0.89, 0.85, 0.80, 0.75, 0]; % Added final 0
            cS.initial_pop = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, ...
                             64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];
            cS.initial_pop = cS.initial_pop / sum(cS.initial_pop);
            cS.bgp_tolerance = 1e-5; cS.bgp_window = 15;
            cS.max_periods = 600; cS.nPeriods = cS.max_periods;
            cS.pop_growth_path = zeros(1, cS.max_periods + 1); % Zero growth default for SS
            cS.alpha = 0.36; cS.delta = 0.08; cS.growth_rate_tech = 0.00;
            cS.theta = 0.04; cS.gamma = 2.0;
            cS.h = [2.2, 2.5, 2.7, 2.9, 3.0, 3.1, 3.1, 3.0, 0, 0, 0, 0, 0, 0, 0, 0];
            cS.tau_g = 0.16; cS.lambda = 0.45;
            cS.pension_return_premium = 0.005; cS.pension_withdrawal_rate = 0.15; cS.pension_withdrawal_tax = 0.03;
            cS.qMin = 0.0; cS.qMax = 0.2;
            cS.rho = 0.96; cS.sigma_e = 0.12; cS.sigmaInitial_e = 0.36; cS.nE = 3; % Use nE=3 for faster debug cycles initially

            % --- Asset Grids (Crucial - USE LARGE VALUES) ---
            cS.nK = 100; cS.kMin = 1e-6;
            cS.kMax = 10000.0; % Keep large value from successful previous test
            cS.kGridV = main_olg_v1_utils.createGrid(cS.kMin, cS.kMax, cS.nK, 'log');
            cS.nP = 100; cS.pMin = 0.0;
            cS.pMax = 10000.0; % Keep large value
            cS.pGridV = main_olg_v1_utils.createGrid(cS.pMin, cS.pMax, cS.nP, 'linear');
            if cS.debugMode == 1, fprintf('  [DEBUG] Grids: nK=%d (Max=%.1f), nP=%d (Max=%.1f), nE=%d\n', cS.nK, cS.kMax, cS.nP, cS.pMax, cS.nE); end
            % --- End Asset Grids ---

            cS.nS_sim = 20000;
            cS.fmincon_options = optimoptions('fmincon', 'Algorithm','sqp', 'Display','none', 'OptimalityTolerance',1e-6, 'StepTolerance',1e-6, 'ConstraintTolerance',1e-7, 'MaxIterations',200, 'MaxFunctionEvaluations',1000);
            cS.interpMethod = 'linear';
            cS.tolLevel_eq = 1e-6; cS.maxIter_eq = 200; cS.dampK = 0.15; cS.dampT = 0.08; cS.dampD = 0.04;
            cS.bequest_transfer_frac = 1.0;
            if cS.debugMode == 1, fprintf('[DEBUG] Parameter Initialization Complete.\n'); end
        end

        % --- createGrid --- (No debug prints needed)
        function grid = createGrid(min_val, max_val, n, type)
             if strcmp(type,'linear'), grid=linspace(min_val,max_val,n); elseif strcmp(type,'log'), if min_val<=0,min_val=1e-9;end; log_min=log(min_val);log_max=log(max_val);log_grid=linspace(log_min,log_max,n);grid=exp(log_grid); grid(1)=min_val; else error('Invalid grid type'); end
         end

        % --- detectSteadyStatePopulation --- (Keep figure for debug)
        function [Z_ss, dependency_ratio, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
             actual_periods=size(popS.Z,2); bgp_reached=false; bgp_period=actual_periods;
             if actual_periods<cS.bgp_window+2
                 if cS.debugMode == 1, fprintf('[DEBUG] Pop sim %d too short for SS check.\n',actual_periods); end
             else
                 psr=popS.Z./sum(popS.Z,1); psc=zeros(actual_periods-1,1); for t=2:actual_periods,psc(t-1)=norm(psr(:,t)-psr(:,t-1));end; drh=zeros(actual_periods,1); for t=1:actual_periods,wp=sum(popS.Z(1:cS.workingAgeMaxIdx,t));rp=sum(popS.Z(cS.retirementAgeIdx:end,t)); if wp>0,drh(t)=rp/wp;else drh(t)=inf;end;end; drc=abs(diff(drh)); ltp=log(popS.totalPop);ltp(~isfinite(ltp))=NaN; grh=diff(ltp);if isempty(grh),grh=NaN;end; grc=abs(diff(grh));if isempty(grc),grc=NaN;end; sst=cS.bgp_window+2;
                 if actual_periods>=sst, for t=sst:actual_periods
                     ss=all(psc(t-cS.bgp_window:t-1)<cS.bgp_tolerance); igs=max(1,t-cS.bgp_window-1);ige=min(length(grc),t-2); gs=false; if igs<=ige, vgc=grc(igs:ige);vgc=vgc(isfinite(vgc));if isempty(vgc),gs=true;else gs=all(vgc<cS.bgp_tolerance);end; else gs=true;end; vdc=drc(t-cS.bgp_window:t-1);vdc=vdc(isfinite(vdc));if isempty(vdc),ds=true;else ds=all(vdc<cS.bgp_tolerance);end;
                     if ss&&gs&&ds, bgp_reached=true;bgp_period=t;break;end; end; end; end
             if bgp_reached, Z_ss=popS.Z(:,bgp_period);dependency_ratio=drh(bgp_period); fprintf('Pop SS reached at period %d\n',bgp_period); else fprintf('Pop SS not reached, using period %d\n',actual_periods); Z_ss=popS.Z(:,actual_periods); if isempty(popS.dependencyRatio), wp=sum(Z_ss(1:cS.workingAgeMaxIdx));rp=sum(Z_ss(cS.retirementAgeIdx:end));if wp>0,dependency_ratio=rp/wp;else dependency_ratio=inf;end; else dependency_ratio=popS.dependencyRatio(actual_periods); end; bgp_period=actual_periods; end
             % Keep figure always, useful visualization
             figure;bar(1:cS.aD,Z_ss/sum(Z_ss)*100);title('SS Pop Structure (%)');xlabel('Age Group');grid on; drawnow;
         end

        % --- initPopulation --- (Keep figure)
        function popS = initPopulation(cS)
            popS.Z=zeros(cS.aD,1);popS.Z(:,1)=cS.initial_pop(:)*100;popS.totalPop=sum(popS.Z(:,1));popS.ageDist=popS.Z(:,1)/popS.totalPop;popS.x=cS.pop_growth_path;
            figure;bar(1:cS.aD,popS.ageDist*100);title('Initial Pop Structure (%)');xlabel('Age Group');grid on; drawnow;
        end

        % --- populationDynamics --- (Keep figures)
        function popS = populationDynamics(popS, cS)
             if length(popS.x)<cS.max_periods+1,warning('popS.x short');popS.x(length(popS.x)+1:cS.max_periods+1)=popS.x(end);end; Z_new=zeros(cS.aD,cS.max_periods);tp_new=zeros(1,cS.max_periods);ad_new=zeros(cS.aD,cS.max_periods); Z_new(:,1)=popS.Z(:,1);tp_new(1)=popS.totalPop(1);ad_new(:,1)=popS.ageDist(:,1); fprintf('Pop dynamics start (Max=%d)...\n',cS.max_periods);bgr=false;ap=cS.max_periods;
             for t=2:cS.max_periods, if mod(t,100)==0||t==2,fprintf('  Pop period %d\n',t);end; gr=popS.x(t);Z_new(1,t)=Z_new(1,t-1)*(1+gr); for a=2:cS.aD,Z_new(a,t)=Z_new(a-1,t-1)*cS.beta(a-1);end; tp_new(t)=sum(Z_new(:,t));if tp_new(t)>0,ad_new(:,t)=Z_new(:,t)/tp_new(t);else ad_new(:,t)=0;end;
                 if t>=cS.bgp_window+2, Zwin=Z_new(:,t-cS.bgp_window:t);tpwin=tp_new(t-cS.bgp_window:t);if any(tpwin<=0),continue;end; psr=Zwin./sum(Zwin,1);psc=zeros(cS.bgp_window,1);for i=1:cS.bgp_window,psc(i)=norm(psr(:,i+1)-psr(:,i));end; drhw=zeros(cS.bgp_window+1,1);for i=0:cS.bgp_window,wp=sum(Zwin(1:cS.workingAgeMaxIdx,i+1));rp=sum(Zwin(cS.retirementAgeIdx:end,i+1));if wp>0,drhw(i+1)=rp/wp;else drhw(i+1)=inf;end;end;drcw=abs(diff(drhw)); ltpw=log(tpwin);ltpw(~isfinite(ltpw))=NaN;grhw=diff(ltpw);grcw=abs(diff(grhw));if isempty(grcw),grcw=NaN;end; ss=all(psc<cS.bgp_tolerance);vgc=grcw(isfinite(grcw));if isempty(vgc),gs=true;else gs=all(vgc<cS.bgp_tolerance);end;vdc=drcw(isfinite(drcw));if isempty(vdc),ds=true;else ds=all(vdc<cS.bgp_tolerance);end;
                 if ss&&gs&&ds,fprintf('\nPop SS reached at period %d\n',t);bgr=true;ap=t;break;end;end;end
             popS.Z=Z_new(:,1:ap);popS.totalPop=tp_new(1:ap);popS.ageDist=ad_new(:,1:ap);popS.x=popS.x(1:ap+1); drh=zeros(ap,1);for th=1:ap,wp=sum(popS.Z(1:cS.workingAgeMaxIdx,th));rp=sum(popS.Z(cS.retirementAgeIdx:end,th));if wp>0,drh(th)=rp/wp;else drh(th)=inf;end;end;popS.dependencyRatio=drh; fprintf('Pop dynamics end. Actual periods: %d\n',ap);if ~bgr,fprintf('Warning: Pop SS not reached.\n');end;
             figure;plot(1:ap,popS.dependencyRatio,'-o');title('Dep Ratio');grid on;drawnow; figure;plot(1:ap,popS.totalPop);title('Total Pop');grid on;drawnow;
         end

        % --- EarningProcess_olgm --- (Conditional debug print)
        function [leLogGridV, leTrProbM, leProb1V] = EarningProcess_olgm(cS)
             m=3;uM=0;uV=cS.sigma_e^2/(1-cS.rho^2);uS=sqrt(uV);zMn=uM-m*uS;zMx=uM+m*uS;leLogGridV=linspace(zMn,zMx,cS.nE);leTrProbM=zeros(cS.nE,cS.nE);gW=(zMx-zMn)/(cS.nE-1);for i=1:cS.nE,for j=1:cS.nE,cMn=cS.rho*leLogGridV(i);if j==1,leTrProbM(i,j)=normcdf((leLogGridV(j)-cMn+gW/2)/cS.sigma_e);elseif j==cS.nE,leTrProbM(i,j)=1-normcdf((leLogGridV(j)-cMn-gW/2)/cS.sigma_e);else leTrProbM(i,j)=normcdf((leLogGridV(j)-cMn+gW/2)/cS.sigma_e)-normcdf((leLogGridV(j)-cMn-gW/2)/cS.sigma_e);end;end;rs=sum(leTrProbM(i,:));if rs>1e-9,leTrProbM(i,:)=leTrProbM(i,:)/rs;else leTrProbM(i,:)=1/cS.nE;end;end;leProb1V=zeros(1,cS.nE);for j=1:cS.nE,if j==1,leProb1V(j)=normcdf((leLogGridV(j)+gW/2)/cS.sigmaInitial_e);elseif j==cS.nE,leProb1V(j)=1-normcdf((leLogGridV(j)-gW/2)/cS.sigmaInitial_e);else leProb1V(j)=normcdf((leLogGridV(j)+gW/2)/cS.sigmaInitial_e)-normcdf((leLogGridV(j)-gW/2)/cS.sigmaInitial_e);end;end;leProb1V=leProb1V/sum(leProb1V);
             if cS.debugMode == 1, fprintf('[DEBUG] Earning process generated.\n'); end
        end

        % --- computeAggLabor --- (Conditional debug print)
        function L = computeAggLabor(popS_Z_ss, cS, paramS)
             L_sum=0;if isempty(paramS.leTrProbM)||size(paramS.leTrProbM,1)~=cS.nE||size(paramS.leTrProbM,2)~=cS.nE,error('Labor matrix invalid.');end;try [evc,eva]=eig(paramS.leTrProbM');eva_d=diag(eva);[~,idx]=min(abs(eva_d-1));pi_e=evc(:,idx);pi_e=real(pi_e);pi_e(pi_e<0)=0;if sum(pi_e)>1e-9,pi_e=pi_e/sum(pi_e);else warning('Using uniform labor shock dist.');pi_e=ones(cS.nE,1)/cS.nE;end;catch ME,warning('Error computing labor shock dist: %s. Using uniform.',ME.message);pi_e=ones(cS.nE,1)/cS.nE;end;mean_raw_e_ss=sum(paramS.leGridV(:).*pi_e(:));
             if cS.debugMode == 1, fprintf('  [DEBUG] computeAggLabor: mean_raw_e_ss = %.4f\n', mean_raw_e_ss); end
             for a=1:cS.workingAgeMaxIdx,L_a=popS_Z_ss(a)*paramS.ageEffV(a)*mean_raw_e_ss;L_sum=L_sum+L_a;end;L=L_sum;
        end

        % --- computePrices --- (Conditional debug prints)
        function [Y, R, w, b, deficit, public_pension_revenue, public_pension_expenditure] = computePrices(K, L, D, Z_ss, cS)
            K=max(K,1e-6);L=max(L,1e-6);Y=K^cS.alpha*L^(1-cS.alpha);
            r_calc = cS.alpha*(Y/K) - cS.delta; w=(1-cS.alpha)*(Y/L); R_unconstrained=1+r_calc; R=max(R_unconstrained,1+1e-6);
            if cS.debugMode == 1
                fprintf('  [DEBUG] computePrices: K=%.2f, L=%.2f, Y=%.2f, K/L=%.2f, Y/K=%.4f\n', K, L, Y, K/L, Y/K);
                fprintf('    Raw r = %.5f (MPK=%.5f), Constrained R = %.5f\n', r_calc, cS.alpha*(Y/K), R);
            end
            wp_c=sum(Z_ss(1:cS.workingAgeMaxIdx));rp_c=sum(Z_ss(cS.retirementAgeIdx:end)); pub_rev=cS.tau_g*w*L;if wp_c>0,avg_w_pwp=(w*L)/wp_c;else avg_w_pwp=0;end; pen_bpr=cS.lambda*avg_w_pwp;pen_bpr=max(0,pen_bpr); pub_exp=pen_bpr*rp_c;pub_exp=max(0,pub_exp); deficit=pub_exp-pub_rev;b=pen_bpr;
            if cS.debugMode == 1, fprintf('    w=%.4f, b=%.4f, deficit=%.4f (Rev=%.2f, Exp=%.2f)\n', w, b, deficit, pub_rev, pub_exp); end
            public_pension_revenue=pub_rev; public_pension_expenditure=pub_exp;
        end

        % =====================================================================
        % ================= solveHouseholdProblem (Conditional DEBUG) =========
        % =====================================================================
        function [cPolM, kPolM, pPolM, qPolM, valueM] = solveHouseholdProblem(R, w, T, b, cS, paramS)
            beta_vfi = 1 / (1 + cS.theta); valueM = zeros(cS.nK, cS.nP, cS.nE, cS.aD); cPolM = valueM; kPolM = valueM; pPolM = valueM; qPolM = valueM; penaltyValue = -1e12;
            if cS.debugMode == 1, fprintf('[DEBUG] 开始求解家户问题... R=%.5f, w=%.4f, T=%.4f, b=%.4f\n', R, w, T, b); end
            nStatesPerAge = cS.nK * cS.nP * cS.nE;

            for a = cS.aD:-1:1
                 startTimeAge = tic;
                 if cS.debugMode == 1, fprintf('  [DEBUG] Solving age %d/%d (%d states)... ', a, cS.aD, nStatesPerAge); end
                 is_working_age = (a <= cS.workingAgeMaxIdx); pension_income_public_age = (1-is_working_age)*b; R_pension_accum = R+cS.pension_return_premium; R_pension_accum=max(R_pension_accum,1.001); current_T_age = (a==1)*T; valueInterp_nextAge=cell(cS.nE,1);
                 if a<cS.aD, next_age=a+1; for e_next=1:cS.nE, Vsn=valueM(:,:,e_next,next_age); fm=isfinite(Vsn); if any(fm(:)), mfv=min(Vsn(fm)); Vsn(~fm)=min(penaltyValue/2,mfv-1e10); else Vsn(:)=penaltyValue/2; end; try valueInterp_nextAge{e_next}=griddedInterpolant({cS.kGridV,cS.pGridV},Vsn,cS.interpMethod,'linear'); catch ME, error('Interp create failed a %d->%d, e %d: %s',a,next_age,e_next,ME.message); end; end; end
                 valueM_a_lin=zeros(nStatesPerAge,1); cPolM_a_lin=zeros(nStatesPerAge,1); kPolM_a_lin=zeros(nStatesPerAge,1); pPolM_a_lin=zeros(nStatesPerAge,1); qPolM_a_lin=zeros(nStatesPerAge,1);
                 kGV=cS.kGridV; pGV=cS.pGridV; leGV=paramS.leGridV; aEV=paramS.ageEffV; fminOpts=cS.fmincon_options;
                 kMin=cS.kMin; kMax=cS.kMax; qMin=cS.qMin; qMax=cS.qMax; pMin=cS.pMin; tauG=cS.tau_g; gam=cS.gamma; pWR=cS.pension_withdrawal_rate; pWT=cS.pension_withdrawal_tax;

                 parfor state_idx = 1:nStatesPerAge
                     e=floor((state_idx-1)/(cS.nK*cS.nP))+1; ki=mod(floor((state_idx-1)/cS.nP),cS.nK)+1; pi=mod(state_idx-1,cS.nP)+1;
                     ck=kGV(ki); cp=pGV(pi); effL=is_working_age*leGV(e)*aEV(a); wie=w*effL;
                     vOpt=penaltyValue; cOpt=1e-6; kpOpt=kMin; ppOpt=pMin; qOpt=qMin;

                     if a==cS.aD, C=(1+R)*ck+cp+pension_income_public_age; if C>1e-6, vOpt=main_olg_v1_utils.utilityFunction(C,gam); cOpt=C; kpOpt=0; ppOpt=0; qOpt=0; end
                     else
                         if is_working_age
                             lb=[kMin,qMin]; ub=[kMax,qMax]; objW=@(x)main_olg_v1_utils.neg_bellman_worker(x,ck,cp,e,a,R,R_pension_accum,w,wie,current_T_age,beta_vfi,cS,paramS,valueInterp_nextAge); x0=[max(kMin,min(kMax,ck)),qMin];
                             try [xf,fv,flg]=fmincon(objW,x0,[],[],[],[],lb,ub,[],fminOpts); catch, xf=[kMin,qMin];fv=-penaltyValue;flg=-99; end
                             if flg>0, kpr=xf(1); qr=xf(2); pcr=qr*wie; ppr=(cp+pcr)*R_pension_accum; ppr=max(pMin,ppr); Res=(1+R)*ck+wie*(1-tauG)+current_T_age; Cr=Res-kpr-pcr;
                                 if Cr>1e-6&&isfinite(fv)&&fv<abs(penaltyValue), vOpt=-fv; cOpt=Cr; kpOpt=kpr; ppOpt=ppr; qOpt=qr;
                                     if cS.debugMode == 1 && state_idx == 1, fprintf('\n  DEBUG HH WORKER (a=%d, kIdx=%d, pIdx=%d, e=%d):\n    k=%.2f, p=%.2f, effL=%.2f, Resources=%.2f\n    fmincon: fval=%.4e, flag=%d\n    Optimal: k''=%.4f, q=%.6f, p''=%.4f, C=%.4f, V=%.4e\n    PensionContrib=%.4f\n', a,ki,pi,e,ck,cp,effL,Res,fv,flg,kpOpt,qOpt,ppOpt,cOpt,vOpt,pcr); end
                                 end; end
                         else % Retirement
                             lb=kMin; ub=kMax; pae=cp*(1-pWR); pae=max(0,pae); ppf=pae; ppf=max(pMin,ppf); pwi=cp*pWR*(1-pWT); Dir=(1+R)*ck+pension_income_public_age+pwi;
                             if Dir-kMin>1e-6, objR=@(x)main_olg_v1_utils.neg_bellman_retiree(x,ck,cp,ppf,e,a,R,pension_income_public_age,Dir,beta_vfi,cS,paramS,valueInterp_nextAge); x0=max(kMin,min(kMax,ck));
                                 try [xf,fv,flg]=fmincon(objR,x0,[],[],[],[],lb,ub,[],fminOpts); catch, xf=kMin;fv=-penaltyValue;flg=-98; end
                                 if flg>0, kpr=xf(1); Cr=Dir-kpr;
                                     if Cr>1e-6&&isfinite(fv)&&fv<abs(penaltyValue), vOpt=-fv; cOpt=Cr; kpOpt=kpr; ppOpt=ppf; qOpt=0;
                                         if cS.debugMode == 1 && state_idx == 1, fprintf('\n  DEBUG HH RETIREE (a=%d, kIdx=%d, pIdx=%d, e=%d):\n    k=%.2f, p=%.2f, DisposableInc=%.2f, PWInc=%.2f\n    fmincon: fval=%.4e, flag=%d\n    Optimal: k''=%.4f, p''(fixed)=%.4f, C=%.4f, V=%.4e\n', a,ki,pi,e,ck,cp,Dir,pwi,fv,flg,kpOpt,ppOpt,cOpt,vOpt); end
                                     end; end
                             else kpOpt=kMin; ppOpt=ppf; qOpt=0; end
                         end % End working age check
                     end % End last period check
                     valueM_a_lin(state_idx)=vOpt; cPolM_a_lin(state_idx)=cOpt; kPolM_a_lin(state_idx)=kpOpt; pPolM_a_lin(state_idx)=ppOpt; qPolM_a_lin(state_idx)=qOpt;
                 end % --- End parfor ---

                 valueM_r=reshape(valueM_a_lin,[cS.nK,cS.nP,cS.nE]); valueM(:,:,:,a)=valueM_r; cPolM_r=reshape(cPolM_a_lin,[cS.nK,cS.nP,cS.nE]); cPolM(:,:,:,a)=cPolM_r;
                 kPolM_r=reshape(kPolM_a_lin,[cS.nK,cS.nP,cS.nE]); kPolM(:,:,:,a)=kPolM_r; pPolM_r=reshape(pPolM_a_lin,[cS.nK,cS.nP,cS.nE]); pPolM(:,:,:,a)=pPolM_r;
                 qPolM_r=reshape(qPolM_a_lin,[cS.nK,cS.nP,cS.nE]); qPolM(:,:,:,a)=qPolM_r;

                 if cS.debugMode == 1
                     if is_working_age, all_q=qPolM_r(:); vq=isfinite(all_q); avq=0; if any(vq),avq=mean(all_q(vq));end; npos=sum(all_q(vq)>1e-4); fprintf('\n  DEBUG HH (Age %d Summary): Avg q=%.6f, States q>1e-4=%d/%d (%.2f%%)\n',a,avq,npos,nStatesPerAge,100*npos/nStatesPerAge);
                     else fprintf('\n  DEBUG HH (Age %d Summary): Retirement Age (q=0)\n', a); end
                 end
                 elapsedTimeAge = toc(startTimeAge); if cS.debugMode == 1, fprintf(' Done in %.2f seconds.\n', elapsedTimeAge); end
            end % --- End age loop ---
            if cS.debugMode == 1, fprintf('[DEBUG] 家户问题求解完成。\n'); end
        end

        % =====================================================================
        % ======================= Helper Static Methods =======================
        % =====================================================================
        % --- neg_bellman_worker --- (No changes, prints handled in main func)
        function neg_value = neg_bellman_worker(x, currentK, currentP, e, a, R, R_pension_accum, w, wage_income_eff, current_T_age, beta_vfi, cS, paramS, valueInterp_nextAge)
             k_prime=x(1); q=x(2); penaltyValue=1e12; pension_contrib=q*wage_income_eff; p_prime=(currentP+pension_contrib)*R_pension_accum; p_prime=max(cS.pMin, p_prime);
             Resources=(1+R)*currentK+wage_income_eff*(1-cS.tau_g)+current_T_age; C=Resources-k_prime-pension_contrib;
             if C<=1e-6, neg_value=penaltyValue+(1e-6-C)*penaltyValue; return; end
             u=main_olg_v1_utils.utilityFunction(C, cS.gamma); if ~isfinite(u), neg_value=penaltyValue; return; end
             expected_future_v=0; next_age=a+1; if next_age<=cS.aD
                  k_interp_pt=max(cS.kGridV(1),min(cS.kGridV(end),k_prime)); p_interp_pt=max(cS.pGridV(1),min(cS.pGridV(end),p_prime));
                  for e_next=1:cS.nE, prob_e_next=paramS.leTrProbM(e,e_next); if prob_e_next>0
                      v_next_interp=-penaltyValue/2; if ~isempty(valueInterp_nextAge{e_next}), try v_next_interp=valueInterp_nextAge{e_next}(k_interp_pt,p_interp_pt); catch, end; end
                      if ~isfinite(v_next_interp), v_next_interp=-penaltyValue/2; end; expected_future_v=expected_future_v+prob_e_next*v_next_interp; end; end; end
             if ~isfinite(expected_future_v), expected_future_v=-penaltyValue/2; end
             survival_prob=cS.beta(a); total_value=u+beta_vfi*survival_prob*expected_future_v;
             if ~isfinite(total_value), neg_value=penaltyValue; else neg_value=-total_value; end
         end

        % --- neg_bellman_retiree --- (No changes, prints handled in main func)
        function neg_value = neg_bellman_retiree(x, currentK, currentP, p_prime_fixed, e, a, R, pension_income_public_age, DisposableIncome_Retire, beta_vfi, cS, paramS, valueInterp_nextAge)
             k_prime=x(1); penaltyValue=1e12; C=DisposableIncome_Retire-k_prime;
             if C<=1e-6, neg_value=penaltyValue+(1e-6-C)*penaltyValue; return; end
             u=main_olg_v1_utils.utilityFunction(C, cS.gamma); if ~isfinite(u), neg_value=penaltyValue; return; end
             expected_future_v=0; next_age=a+1; if next_age<=cS.aD
                  k_interp_pt=max(cS.kGridV(1),min(cS.kGridV(end),k_prime)); p_interp_pt=max(cS.pGridV(1),min(cS.pGridV(end),p_prime_fixed));
                  for e_next=1:cS.nE, prob_e_next=paramS.leTrProbM(e,e_next); if prob_e_next>0
                      v_next_interp=-penaltyValue/2; if ~isempty(valueInterp_nextAge{e_next}), try v_next_interp=valueInterp_nextAge{e_next}(k_interp_pt,p_interp_pt); catch, end; end
                      if ~isfinite(v_next_interp), v_next_interp=-penaltyValue/2; end; expected_future_v=expected_future_v+prob_e_next*v_next_interp; end; end; end
             if ~isfinite(expected_future_v), expected_future_v=-penaltyValue/2; end
             survival_prob=cS.beta(a); total_value=u+beta_vfi*survival_prob*expected_future_v;
             if ~isfinite(total_value), neg_value=penaltyValue; else neg_value=-total_value; end
         end
        % =====================================================================
        % ===================== End Helper Static Methods =====================
        % =====================================================================

        % --- computeAggregates --- (Conditional DEBUG prints and save)
        function [KModel, PModel, TModel, DModel, agg_C, agg_pension_contrib_private, sim_k_val, sim_p_val, sim_e_idx_val] = computeAggregates(kPolM, pPolM, qPolM, cPolM, Z_ss, T_Guess, D_Guess, deficit_flow, R, w, cS, paramS)
             if cS.debugMode==1, fprintf('[DEBUG] 开始计算宏观总量... R=%.5f, w=%.4f, T_G=%.4f, D_G=%.2f\n', R, w, T_Guess, D_Guess); end
             sim_k_val=zeros(cS.nS_sim,cS.aD); sim_p_val=zeros(cS.nS_sim,cS.aD); sim_c_val=zeros(cS.nS_sim,cS.aD); sim_q_val=zeros(cS.nS_sim,cS.aD); sim_e_idx_val=zeros(cS.nS_sim,cS.aD);
             sim_k_val(:,1)=cS.kMin+T_Guess; sim_p_val(:,1)=cS.pMin; cdf1V=cumsum(paramS.leProb1V); randV_e=rand(cS.nS_sim,1); for i=1:cS.nS_sim, sim_e_idx_val(i,1)=sum(randV_e(i)>cdf1V)+1; end

             kPolInterp=cell(cS.nE,cS.aD); pPolInterp=cell(cS.nE,cS.aD); qPolInterp=cell(cS.nE,cS.aD); cPolInterp=cell(cS.nE,cS.aD);
             if cS.debugMode==1, fprintf('  [DEBUG] 创建聚合插值器...\n'); end
             for a=1:cS.aD, for e=1:cS.nE, ks=kPolM(:,:,e,a); ks(~isfinite(ks))=cS.kMin; ps=pPolM(:,:,e,a); ps(~isfinite(ps))=cS.pMin; qs=qPolM(:,:,e,a); qs(~isfinite(qs))=0; cs=cPolM(:,:,e,a); cs(~isfinite(cs)|cs<=0)=1e-6; try kPolInterp{e,a}=griddedInterpolant({cS.kGridV,cS.pGridV},ks,cS.interpMethod,'linear'); pPolInterp{e,a}=griddedInterpolant({cS.kGridV,cS.pGridV},ps,cS.interpMethod,'linear'); qPolInterp{e,a}=griddedInterpolant({cS.kGridV,cS.pGridV},qs,cS.interpMethod,'linear'); cPolInterp{e,a}=griddedInterpolant({cS.kGridV,cS.pGridV},cs,cS.interpMethod,'linear'); catch ME, error('Interp create failed a %d, e %d: %s',a,e,ME.message); end; end; end;
             if cS.debugMode==1, fprintf('  [DEBUG] 插值器创建完毕.\n'); end

             % --- Conditional SAVE pPolM for Debugging ---
             if cS.debugMode == 1
                 debugFolderName = 'tmp'; if ~exist(debugFolderName, 'dir'), mkdir(debugFolderName); fprintf('  [DEBUG] Created folder: %s\n', debugFolderName); end
                 saveFilename = fullfile(debugFolderName, 'policy_debug.mat');
                 try save(saveFilename, 'kPolM', 'pPolM', 'qPolM', 'cPolM'); fprintf('  [DEBUG] Saved policy functions to %s\n', saveFilename);
                 catch ME_save, warning('  [DEBUG] Failed to save policy functions: %s\n', ME_save.message); end
             end
             % --- End Conditional SAVE ---

             if cS.debugMode==1, fprintf('  [DEBUG] 模拟 %d 个体路径...\n', cS.nS_sim); end
             debug_i=1; debug_a_work=2; debug_a_retire=10; % Indices for debug prints inside sim loop

             for a=1:cS.aD, k_int=kPolInterp(:,a); p_int=pPolInterp(:,a); q_int=qPolInterp(:,a); c_int=cPolInterp(:,a); leTrP=paramS.leTrProbM;
                 for i=1:cS.nS_sim
                     ck=sim_k_val(i,a); cp=sim_p_val(i,a); ce=sim_e_idx_val(i,a);
                     if ce<1||ce>cS.nE||~isfinite(ce), sim_c_val(i,a)=1e-6; sim_q_val(i,a)=0; if a<cS.aD, sim_k_val(i,a+1)=cS.kMin; sim_p_val(i,a+1)=cS.pMin; sim_e_idx_val(i,a+1)=1; end; continue; end
                     ki=max(cS.kGridV(1),min(cS.kGridV(end),ck)); pi=max(cS.pGridV(1),min(cS.pGridV(end),cp));
                     try nk=k_int{ce}(ki,pi); np=p_int{ce}(ki,pi); cc=c_int{ce}(ki,pi); cq=q_int{ce}(ki,pi); catch, nk=cS.kMin; np=cS.pMin; cc=1e-6; cq=0; end
                     nk=max(cS.kMin,min(cS.kMax,nk)); if ~isfinite(nk),nk=cS.kMin;end; np=max(cS.pMin,min(cS.pMax,np)); if ~isfinite(np),np=cS.pMin;end; cc=max(1e-6,cc); if ~isfinite(cc),cc=1e-6;end; cq=max(cS.qMin,min(cS.qMax,cq)); if ~isfinite(cq),cq=cS.qMin;end;
                     if cS.debugMode == 1 && i == debug_i && (a == debug_a_work || a == debug_a_retire), fprintf('  DEBUG AGG Sim (i=%d, a=%d): State(k=%.2f, p=%.2f, e=%d) -> Policy(k''=%.2f, p''=%.2f, c=%.2f, q=%.4f)\n', i, a, ck, cp, ce, nk, np, cc, cq); end
                     sim_c_val(i,a)=cc; sim_q_val(i,a)=cq;
                     if a<cS.aD, sim_k_val(i,a+1)=nk; sim_p_val(i,a+1)=np; try cdft=cumsum(leTrP(ce,:)); rn=rand(); sim_e_idx_val(i,a+1)=sum(rn>cdft)+1; catch, sim_e_idx_val(i,a+1)=1; end; if sim_e_idx_val(i,a+1)<1||sim_e_idx_val(i,a+1)>cS.nE||~isfinite(sim_e_idx_val(i,a+1)), sim_e_idx_val(i,a+1)=1; end; end
                 end; end;
             if cS.debugMode==1, fprintf('  [DEBUG] 模拟完毕。\n'); end

             K_tot=0; P_tot=0; C_tot=0; ppc_tot=0; sim_lab=zeros(cS.nS_sim,cS.aD);
             for al=1:cS.workingAgeMaxIdx, for il=1:cS.nS_sim, el=sim_e_idx_val(il,al); if el>=1&&el<=cS.nE, sim_lab(il,al)=paramS.leGridV(el)*paramS.ageEffV(al); else sim_lab(il,al)=0; end; end; end
             total_pop_ss=sum(Z_ss); if total_pop_ss<=0, error('SS pop zero.'); end

             if cS.debugMode == 1, fprintf('  [DEBUG] Aggregating...\n'); end
             for a_agg = 1:cS.aD
                  vk=isfinite(sim_k_val(:,a_agg)); vp=isfinite(sim_p_val(:,a_agg)); vc=isfinite(sim_c_val(:,a_agg))&sim_c_val(:,a_agg)>0;
                  avg_k=0; if any(vk), avg_k=mean(sim_k_val(vk,a_agg)); end; avg_p=0; if any(vp), avg_p=mean(sim_p_val(vp,a_agg)); end; avg_c=0; if any(vc), avg_c=mean(sim_c_val(vc,a_agg)); end
                  if cS.debugMode == 1, fprintf('    DEBUG AGG Age %2d: AvgK=%.2f, AvgP=%.2f, AvgC=%.2f\n', a_agg, avg_k, avg_p, avg_c); end
                  K_tot=K_tot+avg_k*Z_ss(a_agg); P_tot=P_tot+avg_p*Z_ss(a_agg); C_tot=C_tot+avg_c*Z_ss(a_agg);
                  if a_agg <= cS.workingAgeMaxIdx
                      vcidx=isfinite(sim_q_val(:,a_agg))&isfinite(sim_lab(:,a_agg)); avg_q_rate=0; avg_pc_u=0; if any(vcidx), avg_q_rate=mean(sim_q_val(vcidx,a_agg)); avg_pc_u=mean(sim_q_val(vcidx,a_agg).*sim_lab(vcidx,a_agg)); end
                      if cS.debugMode == 1, fprintf('      DEBUG AGG Worker: AvgQ=%.6f, AvgContribUnits=%.6f\n', avg_q_rate, avg_pc_u); end
                      if exist('w','var')&&isfinite(w), ppc_tot=ppc_tot+avg_pc_u*w*Z_ss(a_agg); end
                  elseif a_agg >= cS.retirementAgeIdx && a_agg < cS.aD % Avoid index out of bound for beta
                      avg_p_next_implied = avg_p * (1 - cS.pension_withdrawal_rate);
                      if cS.debugMode == 1, fprintf('      DEBUG AGG Retiree: AvgP=%.2f, Implied Next P=%.2f (using AvgP)\n', avg_p, avg_p_next_implied); end
                  end
             end

             KModel=K_tot+P_tot; PModel=P_tot; agg_C=C_tot; agg_pension_contrib_private=ppc_tot;
             T_tot_beq=0; for ab=1:cS.aD-1, dp=1-cS.beta(ab); vnk=isfinite(sim_k_val(:,ab+1)); vnp=isfinite(sim_p_val(:,ab+1)); vnidx=vnk&vnp; avg_na=0; if any(vnidx), avg_na=mean(sim_k_val(vnidx,ab+1)+sim_p_val(vnidx,ab+1)); end; if isnan(avg_na),avg_na=0;end; T_tot_beq=T_tot_beq+Z_ss(ab)*dp*avg_na; end
             TModel=(T_tot_beq/total_pop_ss)*cS.bequest_transfer_frac;
             DModel=D_Guess*R+deficit_flow;

             if cS.debugMode == 1
                 fprintf('[DEBUG] 宏观总量计算完成。\n');
                 fprintf('  KModel=%.2f, PModel=%.2f, TModel=%.4f, DModel=%.2f, AggC=%.2f, AggPPContrib=%.6f\n', KModel,PModel,TModel,DModel,agg_C, agg_pension_contrib_private);
             end
        end

        % --- computeAverageAssetsByAge --- (No debug prints needed)
        function avg_assets = computeAverageAssetsByAge(sim_k_val, sim_p_val, age, cS)
              if age<1||age>cS.aD,error('Invalid age:%d',age);end;nS=size(sim_k_val,1);if nS<=0,avg_assets=0;return;end;vk=isfinite(sim_k_val(:,age));vp=isfinite(sim_p_val(:,age));avg_k=0;if any(vk),avg_k=mean(sim_k_val(vk,age));end;avg_p=0;if any(vp),avg_p=mean(sim_p_val(vp,age));end;avg_assets=avg_k+avg_p;
        end

        % --- utilityFunction --- (Conditional warning)
        function utility = utilityFunction(c, gamma)
            min_c = 1e-9; utility = -1e12;
            if isscalar(c)
                if c > min_c, if gamma==1, utility=log(c); else utility=(c.^(1-gamma))/(1-gamma); end; if ~isfinite(utility), utility=-1e12; end
                else
                    % if cS.debugMode == 1 % Need cS passed in to make this conditional
                    %     fprintf('Warning: utilityFunction called with C = %.4e <= %.1e\n', c, min_c);
                    % end
                    utility = -1e12 - (min_c-c)*1e12;
                end
            else, error('Utility function expects scalar C'); end
        end

    end % End methods
end % End classdef
% --- END OF FILE main_olg_v1_utils.m ---