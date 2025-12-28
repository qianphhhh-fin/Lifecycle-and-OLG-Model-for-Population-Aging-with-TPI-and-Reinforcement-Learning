import mdptoolbox.example
P, R = mdptoolbox.example.forest()
vi = mdptoolbox.mdp.PolicyIteration(P, R, 0.9)
# mdptoolbox.mdp.FiniteHorizon
vi.run()
print(vi.policy) # result is (0, 0, 0)