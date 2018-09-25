import zmodels.lti as tf
import zmodels.fractional as ne

import numpy as np
import matplotlib.pyplot as plt


# physical constants of the model
K = 5e3
gamma = 0.72
Wu = 4e2
# frequencial bounding of the validity of the fractional derivative
Wb = 1e-3
Wh = 1e7

# fractional derivator
D = ne.oustaloup_approx(gamma,Wb,Wh,25)
w, mag, phase = tf.bode(D)

mag_theo, phase_theo = ne.theoretical_bode(gamma,w)

plt.figure()
plt.semilogx(w, mag)    # Bode magnitude plot
plt.semilogx(w, mag_theo)
plt.grid(True, which="both")
plt.figure()
plt.semilogx(w, phase)  # Bode phase plot
plt.semilogx(w, phase_theo)
plt.grid(True, which="both")

# computation of the impedance model
Z_model = K*(((D/(Wu**gamma))+1)/(D/(Wu**gamma)))

W_z = np.logspace(1,5,base=10.0,num=51)

w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
plt.figure()
plt.semilogx(w_z, mag_z)
plt.grid(True, which="both")

# # show all figures
plt.show()
