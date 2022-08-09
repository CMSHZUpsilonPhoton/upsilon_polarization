import numpy as np
import matplotlib.pyplot as plt

from upsilon_polarization import extreme_scenarios

if __name__ == "__main__":
    N_SAMPLE = 1000000
    cos_angles_uniform = np.random.uniform(-1,+1,N_SAMPLE)
    cos_angles_uniform_weights = np.ones(N_SAMPLE)

    cos_angles_uniform_pol_plus = np.random.uniform(-1,+1,N_SAMPLE)
    cos_angles_uniform_pol_plus_weights = []
    for cos_theta in cos_angles_uniform_pol_plus:
        cos_angles_uniform_pol_plus_weights.append(extreme_scenarios.get_weight(cos_theta, +1))

    cos_angles_uniform_pol_minus = np.random.uniform(-1,+1,N_SAMPLE)
    cos_angles_uniform_pol_minus_weights = []
    for cos_theta in cos_angles_uniform_pol_minus:
        cos_angles_uniform_pol_minus_weights.append(extreme_scenarios.get_weight(cos_theta, -1))


    fig = plt.figure()
    ax = plt.axes()

    N_BINS = 40

    ax.hist(cos_angles_uniform, bins=N_BINS, weights=cos_angles_uniform_weights, label = "Uniform", alpha=1, histtype='step')
    ax.hist(cos_angles_uniform_pol_plus, bins=N_BINS, weights=cos_angles_uniform_pol_plus_weights, label = "Plus", alpha=1, histtype='step')
    ax.hist(cos_angles_uniform_pol_minus, bins=N_BINS, weights=cos_angles_uniform_pol_minus_weights, label = "Minus", alpha=1, histtype='step')

    ax.legend()

    fig.savefig("validation_outputs/extreme_scenarios.pdf")
    fig.savefig("validation_outputs/extreme_scenarios.png")
