import numpy as np
import matplotlib.pyplot as plt
import time
import matlab.engine

from main_olg_v9_utils import OLG_V9_Utils, ModelParameters 

SAVE_FRAC = 0.30

def run_python_naive_simulation(eIdxM_group, M, cS, paramS):
    """[带调试打印功能的版本] 使用幼稚策略在纯Python中运行模拟"""
    print("🐍 2. Running Python naive simulation...")
    nSim = cS.nSim
    aD = cS.aD_new

    kHistM = np.zeros((nSim, aD))
    kPpsHistM = np.zeros((nSim, aD))
    cHistM = np.zeros((nSim, aD))
    cppsHistM = np.zeros((nSim, aD))

    k_next = np.full(nSim, cS.kMin)
    k_pps_next = np.full(nSim, cS.kppsMin)

    pps_return_factor = 1 + ((M['R_k_net_factor'] - 1) + cS.pps_return_rate_premium)

    for a_idx in range(aD): # 0-based index
        k_now = k_next
        k_pps_now = k_pps_next
        kHistM[:, a_idx] = k_now
        kPpsHistM[:, a_idx] = k_pps_now

        k_prime_decisions = np.zeros(nSim)
        c_pps_decisions = np.zeros(nSim)
        c_decisions = np.zeros(nSim)
        
        for i_sim in range(nSim):
            is_debugging = (i_sim == 0 and a_idx < 3)
            
            epsilon_val = paramS.leGridV[eIdxM_group[i_sim, a_idx]]
            
            if is_debugging:
                print(f'\n--- PYTHON DEBUG (i_sim=0, a_idx={a_idx}) ---')
                print(f"State In: k_now={k_now[i_sim]:.8f}, k_pps_now={k_pps_now[i_sim]:.8f}, eps_val={epsilon_val:.8f}")

            # --- Start of expanded naive_policy logic for debugging ---
            c_pps_decision = 0.0
            paramS_hh_step = ModelParameters()
            paramS_hh_step.tau_l = M['tau_l']
            paramS_hh_step.pps_tax_deferral_active = cS.pps_active
            paramS_hh_step.current_pps_withdrawal = 0.0
            if a_idx >= cS.aR_new and cS.pps_active:
                 paramS_hh_step.current_pps_withdrawal = k_pps_now[i_sim] * cS.pps_withdrawal_rate

            # --- 深入 HHIncome_Huggett 进行调试 ---
            # 1. 计算非资本收入
            non_capital_income = 0.0
            actual_pps_expenditure = 0.0
            if a_idx < cS.aR_new:
                age_efficiency = cS.ageEffV_new[a_idx]
                labor_income_gross = M['w_gross'] * age_efficiency * epsilon_val
                labor_income_taxable = labor_income_gross # 幼稚策略不交pps，所以没有税收递延
                labor_income_tax = paramS_hh_step.tau_l * max(0, labor_income_taxable)
                labor_income_net = labor_income_gross - labor_income_tax
                actual_pps_expenditure = c_pps_decision
                non_capital_income = labor_income_net + M['TR_total'] + M['bV_payg'][a_idx]
                if is_debugging:
                    print(f'  [PYTHON HHIncome DEBUG a_idx={a_idx}]')
                    print(f'    Work Period: age_eff={age_efficiency:.6f}, labor_inc_gross={labor_income_gross:.6f}, tax={labor_income_tax:.6f}, net={labor_income_net:.6f}')
            else:
                non_capital_income = M['TR_total'] + M['bV_payg'][a_idx]

            # 2. 计算税后资本收入
            capital_income_net_of_tax = k_now[i_sim] * (M['R_k_net_factor'] - 1)
            if is_debugging:
                print(f'    Capital Income: k_now={k_now[i_sim]:.6f}, R_k_net_factor={M["R_k_net_factor"]:.6f}, cap_inc_net={capital_income_net_of_tax:.6f}')
            
            # 3. 计算总资源
            resources = k_now[i_sim] + capital_income_net_of_tax + non_capital_income - actual_pps_expenditure
            # --- 调试结束 ---

            if is_debugging:
                print(f'HHIncome Out: resources={resources:.8f}, non_cap_inc={non_capital_income:.8f}')

            k_prime_decision = SAVE_FRAC * resources
            c_expend = resources - k_prime_decision
            c_decision = max(cS.cFloor, c_expend / (1 + cS.tau_c))

            if is_debugging:
                print(f'Decision: k_prime={k_prime_decision:.8f}, c_decision={c_decision:.8f}')
                print(f'-------------------------------------------')
            # --- End of expanded logic ---

            k_prime_decisions[i_sim] = k_prime_decision
            c_pps_decisions[i_sim] = c_pps_decision
            c_decisions[i_sim] = c_decision
        
        cHistM[:, a_idx] = c_decisions
        cppsHistM[:, a_idx] = c_pps_decisions

        if a_idx < aD - 1:
            k_next = np.clip(k_prime_decisions, cS.kMin, cS.kMax)
            pps_withdrawal = 0
            if a_idx >= cS.aR_new and cS.pps_active:
                pps_withdrawal = k_pps_now * cS.pps_withdrawal_rate
            k_pps_next_unclamped = (k_pps_now + c_pps_decisions - pps_withdrawal) * pps_return_factor
            k_pps_next = np.clip(k_pps_next_unclamped, cS.kppsMin, cS.kppsMax)
            
    print("🐍 Python simulation finished.")
    return kHistM, kPpsHistM, cHistM, cppsHistM


def run_matlab_naive_simulation(eng, eIdxM_group, M, cS, paramS):
    """通过MATLAB引擎调用幼稚策略模拟器"""
    print("Ⓜ️ 3. Running MATLAB naive simulation...")
    cS_dict = {key: val for key, val in cS.__dict__.items() if not key.startswith('__')}
    paramS_dict = {key: val for key, val in paramS.__dict__.items() if not key.startswith('__')}
    eIdxM_ml = matlab.double((eIdxM_group + 1).tolist()) 
    bV_payg_ml = matlab.double(M['bV_payg'].tolist())

    k_ml, kpps_ml, c_ml, cpps_ml = eng.HHSimulation_olgm_naive(
        eIdxM_ml, M['R_k_net_factor'], M['w_gross'], M['TR_total'],
        bV_payg_ml, paramS_dict, cS_dict, nargout=4)

    print("Ⓜ️ MATLAB simulation finished.")
    return np.array(k_ml), np.array(kpps_ml), np.array(c_ml), np.array(cpps_ml)


def main():
    print("--- Naive Policy Consistency Test (Detailed Debug) ---")
    
    print("🚀 1. Initializing parameters and common inputs...")
    cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
    cS.nSim = 100
    cS.nk = 20
    cS.nkpps = 20
    
    paramS = ModelParameters()
    leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)
    paramS.leGridV = np.exp(leLogGridV)
    paramS.leTrProbM = leTrProbM
    paramS.leProb1V = leProb1V

    eIdxM_group = OLG_V9_Utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, paramS.leProb1V, paramS.leTrProbM)

    M_fixed = {
        'R_k_net_factor': 1.03,
        'w_gross': 1.5,
        'TR_total': 0.0,
        'bV_payg': np.concatenate([np.zeros(cS.aR_new), np.full(cS.aD_new - cS.aR_new, 0.2)]),
        'tau_l': 0.25
    }
    
    paramS.tau_l = M_fixed['tau_l']
    paramS.pps_tax_deferral_active = cS.pps_active 
    
    eng = None
    try:
        print("Starting MATLAB engine...")
        eng = matlab.engine.start_matlab()
        eng.addpath(eng.genpath('.'), nargout=0)

        k_py, kpps_py, c_py, cpps_py = run_python_naive_simulation(eIdxM_group, M_fixed, cS, paramS)
        k_ml, kpps_ml, c_ml, cpps_ml = run_matlab_naive_simulation(eng, eIdxM_group, M_fixed, cS, paramS)

        print("\n📊 4. Comparing results...")
        diff_k = np.abs(k_py - k_ml)
        diff_c = np.abs(c_py - c_ml)
        print(f"Max difference in k:    {np.max(diff_k):.2e}")
        print(f"Max difference in c:    {np.max(diff_c):.2e}")
        
        if np.max(diff_k) < 1e-7 and np.max(diff_c) < 1e-7:
            print("\n✅ SUCCESS: The results are numerically identical!")
        else:
            print("\n❌ FAILED: The results are different. Check DEBUG output.")

        plt.figure(figsize=(14, 8))
        plt.subplot(2, 2, 1)
        plt.plot(np.mean(k_py, axis=0), 'b-', label='Python Sim')
        plt.plot(np.mean(k_ml, axis=0), 'r--', label='MATLAB Sim')
        plt.title('Average Capital Path (k)')
        plt.legend(); plt.grid(True)
        
        plt.subplot(2, 2, 2)
        plt.plot(np.mean(c_py, axis=0), 'b-', label='Python Sim')
        plt.plot(np.mean(c_ml, axis=0), 'r--', label='MATLAB Sim')
        plt.title('Average Consumption Path (c)')
        plt.legend(); plt.grid(True)
        
        plt.subplot(2, 2, 3)
        im = plt.imshow(diff_k, aspect='auto', cmap='viridis', origin='lower'); plt.colorbar(im)
        plt.title('Absolute Difference in Capital (k)')
        
        plt.subplot(2, 2, 4)
        im = plt.imshow(diff_c, aspect='auto', cmap='viridis', origin='lower'); plt.colorbar(im)
        plt.title('Absolute Difference in Consumption (c)')
        
        plt.tight_layout()
        plt.suptitle("Python vs MATLAB Simulation Consistency Check", fontsize=16)
        plt.subplots_adjust(top=0.92)
        plt.show()

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if eng:
            print("Quitting MATLAB engine.")
            eng.quit()

if __name__ == '__main__':
    main()