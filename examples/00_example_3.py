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
print 'found for model %s the parameters: %s'%(experiment_1.model[0],str(experiment_1.param[0]))
# Initial guess for the RCPE model
gamma_0=0.7;
K_0=1e3;
wb_0=1e2*2*math.pi;
experiment_1.least_square_fit('Serial_RCPE', p0=(gamma_0, K_0, wb_0))
print 'found for model %s the parameters: %s'%(experiment_1.model[1],str(experiment_1.param[1]))
# Initial guess for the RCR model
K_0 = 5e4
wn_0 = 1e3
wd_0 = 1e1
experiment_1.least_square_fit('RCR_Non_fractional', p0=(K_0, wn_0, wd_0),bounds=([K_0/100, wn_0/100, wd_0/10],[K_0*100, wn_0*100, wd_0*10]))
print 'found for model %s the parameters: %s'%(experiment_1.model[2],str(experiment_1.param[2]))
# EXPERIMENTAL initial guess for the QQR model
gamma_1 = 0.75
wb_1 = 1e7
gamma_2 = 0.1
wb_2 = 1e9
R = 1e3
D_1 = ne.oustaloup_approx(gamma_1,mod.Wl,mod.Wh,mod.N_zp)
D_2 = ne.oustaloup_approx(gamma_2,mod.Wl,mod.Wh,mod.N_zp)
Z_model = 1/(D_1/(wb_1**gamma_1))+1/(D_2/(wb_2**gamma_2))+R
w_z, mag_z, phase_z = tf.bode(Z_model,w=experiment_1.freq)
experiment_1.least_square_fit('Serial_QQR', p0=(gamma_1, wb_1,gamma_2, wb_2, R),bounds=([0.7, wb_1/1000,0.05, wb_2/1000, R/100],[0.8, wb_1*1000,0.15, wb_2*1000, R*100]))
print 'found for model %s the parameters: %s'%(experiment_1.model[3],str(experiment_1.param[3]))

# ploting the experimental values and results of identification
plt.figure()
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB, label='measurement')
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB_ident[0], label=experiment_1.model[0])
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB_ident[1], label=experiment_1.model[1])
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB_ident[2], label=experiment_1.model[2])
#plt.semilogx(experiment_1.freq, mag_z, label='test QQR')
plt.semilogx(experiment_1.freq, experiment_1.zimp_dB_ident[3], label=experiment_1.model[3])
plt.grid(True, which="both")
plt.legend()
plt.xlabel('frequency (Hz)')
plt.ylabel('Impedance magnitude (dB)')
### show figures
plt.show()