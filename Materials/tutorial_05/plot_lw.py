import json
import matplotlib.pyplot as plt
import numpy as np
ff = './MoS2_lambda.json'
with open(ff,'r') as f:
    q_pts = json.load(f)

first_q = q_pts["1"]

print("Fraction coordinate of the q point:", first_q["xq_frac"])
print("Electronic temperaturs used:", first_q["T"])
print("Results for the first mode:", first_q["1"])
print("Frequency of the second mode in meV:", first_q["2"]["freq"])
for q in list(q_pts.values())[:]:
    for mod in range(1,10):
        plt.plot( q["T"],q[f'{mod}']["gamma_allen"],label=r"$\omega_"+f"{mod}$={q[f'{mod}']['freq']:.0f} (meV)",linestyle='--' )
    plt.xlabel('T (eV)')
    plt.ylabel(r'$\gamma(T)$ (meV)')
    plt.legend(title=f"q={ np.round(np.array(q['xq_frac']),2) }")
    plt.show()
