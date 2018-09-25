import zmodels.MeasSignal as expe
import zmodels.Models as mod
import zmodels.lti as tf
import zmodels.fractional as ne

import matplotlib.pyplot as plt
import math

# loading experimental values
experiment_1 = expe.MeasSignal('test2.csv','freq')


# identifications
experiment_1.identif('Serial_RC','LeastSquare')
K_0 = 1e3
wc_0 = 2e4
test = experiment_1.identif('Serial_RC','MinimizeError',p0 = (K_0,wc_0))

# ploting the experimental values and results of identification
plt.figure()
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB, label='measurement')
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB_ident[0], label='test LS')
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB_ident[1], label='test min')
plt.grid(True, which="both")
plt.legend()
plt.xlabel('frequency (Hz)')
plt.ylabel('Impedance magnitude (dB)')
### show figures
plt.show()