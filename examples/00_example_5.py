import zmodels.lti as tf
l
import numpy as np

def trN(H,n):
	# this function returns the time o
	t, step = tf.step2(H,N=500)
	s_inf = step[-1]
	index5 = np.nonzero( abs(step-s_inf) > (n/100.)*s_inf)
	trn = t[np.max(index5)]
	return trn

# curve on Nb_points points
Nb_points = 1000

# second order parameters
K = 1.0
omega_0 = 1.0
m = np.logspace(-2,2,num=Nb_points)
tps_reduit = np.zeros(Nb_points)

# loop on m parameters for the lti
for i in xrange(len(m)):
	H = tf.lti([K],[1/(omega_0*omega_0),(2*m[i])/(float(omega_0)),1])
	tps_reduit[i] = trN(H,5.0)*omega_0

plt.figure()
plt.loglog(m,tps_reduit)
plt.grid()


# write results in a csv file
file = open('data_tps_reduit_vs_m.csv','w')
file.write('tps_reduit,m\n')
for f in xrange(len(tps_reduit)):
	strcsv = str(tps_reduit[f])+','+str(m[f])+'\n'
	file.write(strcsv)
file.close()

plt.show()