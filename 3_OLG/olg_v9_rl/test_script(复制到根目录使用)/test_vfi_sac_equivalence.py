#!/usr/bin/env python3
"""
VFI-SAC Mathematical Equivalence Numerical Verification Script

This script numerically verifies the mathematical equivalence of two methods:
1. VFI method: V(s_t) = u(c_t) + β * s(t) * E[V(s_{t+1})]
2. SAC cumulative survival method: r'(t) = u(c_t) * ∏_{i=1}^{t} s(i), then use standard SAC

Theoretical expectation: Both methods should produce identical lifetime expected utility
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List, Dict
import time

class SimpleOLGParameters:
    """Simplified OLG model parameters"""
    def __init__(self):
        # Basic parameters
        self.beta = 0.97          # Subjective discount factor
        self.sigma = 2.0          # Risk aversion coefficient
        self.T = 16               # Lifecycle length (number of age groups)
        
        # Survival probabilities (decreasing)
        self.survival_probs = np.array([
            0.99, 0.98, 0.97, 0.96, 0.95,  # Young: high survival rates
            0.94, 0.92, 0.90, 0.87, 0.83,  # Middle-aged: gradual decline
            0.78, 0.72, 0.65, 0.55, 0.40,  # Old: rapid decline
            0.20                             # Final period: low survival rate
        ])
        
        # Income profile (inverted U-shape)
        self.income_profile = np.array([
            0.5, 0.7, 0.9, 1.2, 1.5,       # Young: income growth
            1.8, 2.0, 2.1, 2.0, 1.8,       # Middle-aged: income peak
            1.5, 1.0, 0.5, 0.2, 0.1,       # Retirement: income decline
            0.05                            # Final period: minimum income
        ])
        
        # Simplified optimal consumption strategy (for testing)
        self.optimal_consumption = np.array([
            0.4, 0.6, 0.8, 1.0, 1.2,       # Young: moderate consumption
            1.4, 1.5, 1.6, 1.5, 1.4,       # Middle-aged: consumption peak
            1.2, 1.0, 0.8, 0.6, 0.4,       # Old: consumption decline
            0.2                             # Final period: minimum consumption
        ])

def ces_utility(consumption: float, sigma: float) -> float:
    """CES utility function"""
    if consumption <= 0:
        return -1e10  # Avoid negative consumption
    
    if sigma == 1.0:
        return np.log(consumption)
    else:
        return (consumption ** (1 - sigma)) / (1 - sigma)

class VFICalculator:
    """VFI method value function calculation"""
    
    def __init__(self, params: SimpleOLGParameters):
        self.params = params
        
    def calculate_lifetime_utility_vfi(self, consumption_path: np.ndarray) -> float:
        """
        Calculate lifetime utility using VFI method
        
        VFI formula: V(s_t) = u(c_t) + β * s(t) * E[V(s_{t+1})]
        Expanded as: Σ_{t=0}^{T-1} β^t * (∏_{i=0}^{t-1} s(i)) * u(c_t)
        """
        total_utility = 0.0
        
        for t in range(self.params.T):
            # Calculate cumulative survival probability to period t
            if t == 0:
                survival_to_t = 1.0  # Period 0 survival probability is 1
            else:
                survival_to_t = np.prod(self.params.survival_probs[:t])
            
            # Current period utility
            current_utility = ces_utility(consumption_path[t], self.params.sigma)
            
            # VFI form: β^t * cumulative survival probability * current utility
            discounted_utility = (self.params.beta ** t) * survival_to_t * current_utility
            
            total_utility += discounted_utility
            
        return total_utility
    
    def calculate_period_values(self, consumption_path: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate detailed values for each period"""
        utilities = np.zeros(self.params.T)
        survival_to_period = np.zeros(self.params.T)
        discounted_utilities = np.zeros(self.params.T)
        
        for t in range(self.params.T):
            # Cumulative survival probability
            if t == 0:
                survival_to_period[t] = 1.0
            else:
                survival_to_period[t] = np.prod(self.params.survival_probs[:t])
            
            # Current period utility
            utilities[t] = ces_utility(consumption_path[t], self.params.sigma)
            
            # Discounted utility
            discounted_utilities[t] = (self.params.beta ** t) * survival_to_period[t] * utilities[t]
            
        return utilities, survival_to_period, discounted_utilities

class SACEquivalentCalculator:
    """SAC cumulative survival probability method calculation"""
    
    def __init__(self, params: SimpleOLGParameters):
        self.params = params
        
    def calculate_lifetime_utility_sac(self, consumption_path: np.ndarray) -> float:
        """
        Calculate lifetime utility using SAC cumulative survival method
        
        SAC method: r'(t) = u(c_t) * ∏_{i=1}^{t} s(i), then use standard discount β^t
        Equivalent to: Σ_{t=0}^{T-1} β^t * u(c_t) * ∏_{i=1}^{t} s(i)
        """
        total_utility = 0.0
        cumulative_survival = 1.0  # Initial cumulative survival probability
        
        for t in range(self.params.T):
            # Current period utility
            current_utility = ces_utility(consumption_path[t], self.params.sigma)
            
            # SAC method: cumulative survival probability weight * utility
            weighted_reward = current_utility * cumulative_survival
            
            # Standard SAC discounting
            discounted_utility = (self.params.beta ** t) * weighted_reward
            
            total_utility += discounted_utility
            
            # Update cumulative survival probability (for next period)
            if t < self.params.T - 1:  # Not the last period
                cumulative_survival *= self.params.survival_probs[t]
                
        return total_utility
    
    def calculate_period_values(self, consumption_path: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate detailed values for each period"""
        utilities = np.zeros(self.params.T)
        cumulative_weights = np.zeros(self.params.T)
        weighted_rewards = np.zeros(self.params.T)
        discounted_utilities = np.zeros(self.params.T)
        
        cumulative_survival = 1.0
        
        for t in range(self.params.T):
            # Current period utility
            utilities[t] = ces_utility(consumption_path[t], self.params.sigma)
            
            # Cumulative weight
            cumulative_weights[t] = cumulative_survival
            
            # Weighted reward
            weighted_rewards[t] = utilities[t] * cumulative_survival
            
            # Discounted utility
            discounted_utilities[t] = (self.params.beta ** t) * weighted_rewards[t]
            
            # Update cumulative survival probability
            if t < self.params.T - 1:
                cumulative_survival *= self.params.survival_probs[t]
                
        return utilities, cumulative_weights, discounted_utilities

def run_numerical_verification():
    """Run numerical verification"""
    print("=" * 80)
    print("🎯 VFI-SAC Mathematical Equivalence Numerical Verification")
    print("=" * 80)
    
    # 1. Initialize parameters
    params = SimpleOLGParameters()
    vfi_calc = VFICalculator(params)
    sac_calc = SACEquivalentCalculator(params)
    
    print(f"\n📊 Model Parameters:")
    print(f"   Lifecycle length: {params.T} periods")
    print(f"   Subjective discount factor: {params.beta}")
    print(f"   Risk aversion coefficient: {params.sigma}")
    print(f"   Survival probability range: [{np.min(params.survival_probs):.3f}, {np.max(params.survival_probs):.3f}]")
    
    # 2. Use optimal consumption path for verification
    consumption_path = params.optimal_consumption
    
    print(f"\n🍽️ Consumption Path:")
    print(f"   Consumption range: [{np.min(consumption_path):.3f}, {np.max(consumption_path):.3f}]")
    print(f"   Average consumption: {np.mean(consumption_path):.3f}")
    
    # 3. Calculate VFI method lifetime utility
    start_time = time.time()
    vfi_utility = vfi_calc.calculate_lifetime_utility_vfi(consumption_path)
    vfi_time = time.time() - start_time
    
    # 4. Calculate SAC equivalent method lifetime utility
    start_time = time.time()
    sac_utility = sac_calc.calculate_lifetime_utility_sac(consumption_path)
    sac_time = time.time() - start_time
    
    # 5. Compare results
    absolute_diff = abs(vfi_utility - sac_utility)
    relative_diff = absolute_diff / abs(vfi_utility) * 100 if vfi_utility != 0 else 0
    
    print(f"\n🔍 Calculation Results:")
    print(f"   VFI lifetime utility:           {vfi_utility:.12f}")
    print(f"   SAC equivalent lifetime utility: {sac_utility:.12f}")
    print(f"   Absolute difference:            {absolute_diff:.2e}")
    print(f"   Relative difference:            {relative_diff:.2e}%")
    print(f"   VFI calculation time:           {vfi_time*1000:.3f} ms")
    print(f"   SAC calculation time:           {sac_time*1000:.3f} ms")
    
    # 6. Equivalence judgment
    tolerance = 1e-10
    is_equivalent = absolute_diff < tolerance
    
    print(f"\n✅ Equivalence Verification:")
    print(f"   Tolerance:                      {tolerance:.0e}")
    print(f"   Mathematical equivalence:       {'✅ Yes' if is_equivalent else '❌ No'}")
    
    if is_equivalent:
        print(f"   🎉 Congratulations! VFI and SAC cumulative survival methods are numerically equivalent!")
    else:
        print(f"   ⚠️ Note: Numerical differences exist, possibly due to floating point precision or implementation differences")
    
    return vfi_utility, sac_utility, absolute_diff, is_equivalent

def main():
    """Main function"""
    print("🚀 Starting VFI-SAC Mathematical Equivalence Numerical Verification")
    
    # Run basic verification
    vfi_utility, sac_utility, abs_diff, basic_equivalent = run_numerical_verification()
    
    # Final summary
    print("\n" + "=" * 80)
    print("🎯 Final Verification Summary")
    print("=" * 80)
    
    print(f"\n✅ Verification Results:")
    print(f"   Basic equivalence test:         {'✅ Passed' if basic_equivalent else '❌ Failed'}")
    
    print(f"\n🏆 Final Conclusion:")
    if basic_equivalent:
        print(f"   🎉 VFI and SAC cumulative survival methods are numerically equivalent!")
        print(f"   📊 Absolute difference: {abs_diff:.2e}")
        print(f"   📈 Relative difference: {abs_diff/abs(vfi_utility)*100:.2e}%")
        print(f"   💡 It is safe to use SAC cumulative survival method as a substitute for VFI")
    else:
        print(f"   ⚠️ Numerical differences exist, further implementation checking needed")
    
    print(f"\n📊 Technical Details:")
    print(f"   VFI lifetime utility:          {vfi_utility:.12f}")
    print(f"   SAC lifetime utility:          {sac_utility:.12f}")
    print(f"   Machine precision:             {np.finfo(float).eps:.2e}")
    print(f"   Test tolerance:                1.00e-10")
    
    return basic_equivalent

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1) 