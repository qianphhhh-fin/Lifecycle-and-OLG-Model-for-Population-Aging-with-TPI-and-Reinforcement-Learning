import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Tuple, List, Dict

# =============================================================================
#  1. Income Process Simulator (Self-contained)
# =============================================================================
class IncomeProcessSimulator:
    """A dedicated class to simulate the income process from Cocco (2005)."""
    def __init__(self):
        self.tb, self.tr, self.td = 22, 61, 100
        self.aa = -2.170042 + 2.700381
        self.b1, self.b2, self.b3 = 0.16818, -0.0323371 / 10, 0.0019704 / 100
        self.smay, self.smav = np.sqrt(0.169993), np.sqrt(0.112572)
        self.ret_fac, self.tau_y = 0.6827, 0.06
        self.f_y = {age: np.exp(self.aa + self.b1*age + self.b2*age**2 + self.b3*age**3) for age in range(self.tb, self.tr)}

    def simulate_lifecycles(self, num_simulations: int, c_prop: float, alpha: float) -> List[float]:
        """Simulates lifecycles and returns the history of consumption at all timesteps."""
        consumption_history = []
        print(f"Starting simulation of {num_simulations} lifecycles...")
        for _ in range(num_simulations):
            P = np.exp(np.random.normal(0, self.smav))
            P_retirement, W = 0.0, 0.0
            for age in range(self.tb, self.td + 1):
                if age < self.tr:
                    u_shock = np.random.normal(0, self.smay)
                    Y_gross = self.f_y[age] * P * np.exp(u_shock)
                    Y = Y_gross * (1 - self.tau_y)
                    z_shock = np.random.normal(0, self.smav)
                    P *= np.exp(z_shock)
                    if age == self.tr - 1: P_retirement = P
                else:
                    Y = self.ret_fac * P_retirement
                
                X = W + Y
                C = c_prop * X
                consumption_history.append(C)
                S = X - C
                # Simplify: evolve wealth with expected return for a more stable consumption distribution
                portfolio_return = alpha * (1.02 + 0.04) + (1 - alpha) * 1.02
                W = S * portfolio_return
        print("Simulation complete.")
        return consumption_history

# =============================================================================
#  2. Utility and Reward Transformation Functions (Modular)
# =============================================================================
def raw_utility_crra(C: np.ndarray, gamma: float) -> np.ndarray:
    """Calculates raw CRRA utility."""
    C = np.maximum(C, 1e-9) # Avoid numerical errors
    if gamma == 1:
        return np.log(C)
    else:
        return (C**(1 - gamma)) / (1 - gamma)

def transform_reward_affine(U: np.ndarray, a: float, b: float) -> np.ndarray:
    """Applies affine transformation."""
    return U * a + b

def transform_reward_tanh(U: np.ndarray, u_ref: float, sensitivity: float, amplitude: float) -> np.ndarray:
    """Applies tanh transformation."""
    return amplitude * np.tanh((U - u_ref) / sensitivity)

# =============================================================================
#  3. Main Execution Block
# =============================================================================
if __name__ == "__main__":
    
    # --- I. Core Configuration Parameters ---
    GAMMA = 3.84
    NUM_SIMULATIONS = 2000 # Number of lifecycles to simulate
    # Use a relatively conservative consumption proportion for testing to better capture the low-consumption region
    NAIVE_C_PROP = 0.2
    NAIVE_ALPHA = 0.6
    
    # --- II. Simulate and Generate Base Data ---
    simulator = IncomeProcessSimulator()
    consumption_data = np.array(simulator.simulate_lifecycles(NUM_SIMULATIONS, NAIVE_C_PROP, NAIVE_ALPHA))
    utility_data = raw_utility_crra(consumption_data, GAMMA)

    # --- III. Parameter Calibration ---
    print("\n" + "="*50)
    print("Parameter Calibration")
    print("="*50)
    
    # a) Affine Transformation Parameter Calibration
    # -----------------------------------------------
    print("\n--- 1. Affine Transformation ---")
    MIN_C_PERCENTILE = 5     # Use the 5th percentile of consumption distribution as the lower anchor
    REF_C_PERCENTILE = 50     # Use the median of consumption distribution as the reference anchor
    REWARD_FLOOR = -5.0       # Desired reward floor
    REFERENCE_REWARD = 1.0    # Desired reference reward
    
    c_min = np.percentile(consumption_data, MIN_C_PERCENTILE)
    c_ref = np.percentile(consumption_data, REF_C_PERCENTILE)
    
    u_min = raw_utility_crra(c_min, GAMMA)
    u_ref = raw_utility_crra(c_ref, GAMMA)
    
    # Solve for a and b
    affine_a = (REFERENCE_REWARD - REWARD_FLOOR) / (u_ref - u_min)
    affine_b = REWARD_FLOOR - affine_a * u_min
    
    print(f"Consumption Lower Anchor ({MIN_C_PERCENTILE}th percentile): C_min = {c_min:.4f} -> U_min = {u_min:.4f}")
    print(f"Consumption Reference Anchor ({REF_C_PERCENTILE}th percentile): C_ref = {c_ref:.4f} -> U_ref = {u_ref:.4f}")
    print(f"Calculated Affine Transformation Parameters:")
    print(f"  a (scale) = {affine_a:.4f}")
    print(f"  b (shift) = {affine_b:.4f}")

    # b) Tanh Transformation Parameter Calibration
    # --------------------------------------------
    print("\n--- 2. Tanh Transformation ---")
    TANH_AMPLITUDE = 1.5
    # A reasonable heuristic for sensitivity is the difference between reference and min utility
    TANH_SENSITIVITY = (u_ref - u_min)
    
    tanh_u_ref = u_ref # Use median utility as the centering point
    
    print(f"Tanh Transformation Parameters:")
    print(f"  Amplitude = {TANH_AMPLITUDE:.4f}")
    print(f"  U_reference = {tanh_u_ref:.4f}")
    print(f"  Sensitivity = {TANH_SENSITIVITY:.4f}")

    # --- IV. Apply Transformations and Analyze Results ---
    affine_rewards = transform_reward_affine(utility_data, affine_a, affine_b)
    tanh_rewards = transform_reward_tanh(utility_data, tanh_u_ref, TANH_SENSITIVITY, TANH_AMPLITUDE)

    print("\n" + "="*50)
    print("Transformed Reward Statistics")
    print("="*50)
    for name, data in [("Raw Utility", utility_data), ("Affine Transformed Reward", affine_rewards), ("Tanh Transformed Reward", tanh_rewards)]:
        print(f"\n--- {name} ---")
        print(f"  Mean: {np.mean(data):.3f}, Std Dev: {np.std(data):.3f}")
        print(f"  Range: [{np.min(data):.3f}, {np.max(data):.3f}]")

    # --- V. Visualization ---
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, 16))
    gs = fig.add_gridspec(3, 2)
    fig.suptitle('Reward Transformation Parameter Calibration and Effect Analysis', fontsize=22, fontweight='bold')

    # a) Consumption vs. Reward Curves
    ax1 = fig.add_subplot(gs[0, :])
    c_range = np.linspace(np.min(consumption_data), np.percentile(consumption_data, 95), 500)
    u_range = raw_utility_crra(c_range, GAMMA)
    ax1.plot(c_range, u_range, label='Raw CRRA Utility U(C)', color='black', linestyle='--', linewidth=2)
    ax1.plot(c_range, transform_reward_affine(u_range, affine_a, affine_b), label='Affine Transformed Reward R(C)', color='blue', linewidth=3)
    ax1.plot(c_range, transform_reward_tanh(u_range, tanh_u_ref, TANH_SENSITIVITY, TANH_AMPLITUDE), label='Tanh Transformed Reward R(C)', color='red', linewidth=3)
    ax1.scatter([c_min, c_ref], [REWARD_FLOOR, REFERENCE_REWARD], color='blue', s=150, zorder=5, label='Affine Anchor Points')
    ax1.set_title('1. Comparison of Transformation Function Shapes (Consumption vs. Reward)', fontsize=16)
    ax1.set_xlabel('Consumption (C)', fontsize=14)
    ax1.set_ylabel('Reward / Utility Value', fontsize=14)
    ax1.legend(fontsize=12)
    ax1.grid(True, which='both', linestyle=':')
    ax1.set_ylim(REWARD_FLOOR - 1, REFERENCE_REWARD + 2)

    # b) Reward Distributions
    ax2 = fig.add_subplot(gs[1, 0])
    sns.histplot(utility_data, ax=ax2, kde=True, bins=50, color='grey')
    ax2.set_title('2a. Raw Utility Distribution', fontsize=16)
    ax2.set_xlabel('Utility Value')
    
    ax3 = fig.add_subplot(gs[1, 1], sharey=ax2)
    sns.histplot(affine_rewards, ax=ax3, kde=True, bins=50, color='blue', label='Affine Reward')
    sns.histplot(tanh_rewards, ax=ax3, kde=True, bins=50, color='red', label='Tanh Reward')
    ax3.set_title('2b. Transformed Reward Distribution', fontsize=16)
    ax3.set_xlabel('Reward Value')
    ax3.legend()

    # c) Reward Gradient (Sensitivity) dR/dC
    ax4 = fig.add_subplot(gs[2, :])
    delta_c = 1e-4
    u_prime = (raw_utility_crra(c_range + delta_c, GAMMA) - u_range) / delta_c
    affine_grad = affine_a * u_prime
    tanh_grad_numerator = (transform_reward_tanh(u_range + u_prime * delta_c, tanh_u_ref, TANH_SENSITIVITY, TANH_AMPLITUDE) - transform_reward_tanh(u_range, tanh_u_ref, TANH_SENSITIVITY, TANH_AMPLITUDE))
    tanh_grad = tanh_grad_numerator / delta_c
    
    ax4.plot(c_range, affine_grad, label='Affine Transformation Gradient dR/dC', color='blue', linewidth=3)
    ax4.plot(c_range, tanh_grad, label='Tanh Transformation Gradient dR/dC', color='red', linewidth=3)
    ax4.set_title('3. Reward Sensitivity to Consumption (Gradient)', fontsize=16)
    ax4.set_xlabel('Consumption (C)', fontsize=14)
    ax4.set_ylabel('Gradient dR/dC', fontsize=14)
    c_25 = np.percentile(consumption_data, 25)
    c_75 = np.percentile(consumption_data, 75)
    ax4.axvspan(c_25, c_75, color='gray', alpha=0.2, label='Core Consumption Distribution (25%-75%)')
    ax4.legend(fontsize=12)
    ax4.set_yscale('log')
    ax4.set_ylim(bottom=1e-5)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()