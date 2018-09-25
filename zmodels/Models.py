#################################################################################
# Copyright (C) 2016 Florian Kolbl												#
#																				#
# This file is part of ZModels													#
#																				#
# ZModels is free software: you can redistribute it and/or modify				#
# it under the terms of the GNU Lesser General Public License as published by	#
# the Free Software Foundation, either version 3 of the License, or 			#
# (at your option) any later version.											#
#																				#
# ZModels is distributed in the hope that it will be useful,					#
# but WITHOUT ANY WARRANTY; without even the implied warranty of				#
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 					#
# GNU Lesser General Public License for more details.							#
#																				#
# You should have received a copy of the GNU Lesser General Public License 		#
# along with ZModels. If not, see <http://www.gnu.org/licenses/>. 				#
#																				#
# First added:  2016-11-04														#
# Last changed: 2017-12-12														#
#################################################################################
import zmodels.lti as tf
import zmodels.fractional as ne
import zmodels.norm as nrm

# Frequencial bounding of the validity of the fractional derivative
Wl = 1e-3
Wh = 1e7
N_zp = 25

def RC_Zero_Pole_freq(W_z, K, wc, wb):
	'''
	transfer function:
	    1+(s/wb)	1
	 Z=K-------- --------
	     (s/wb)	 1+(s/wc)				'''
	Z_model = K*tf.lti([1/wb,1], [1/wb,0])*tf.lti([1], [1/wc,1])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def RC_Zero_Pole_freq2(W_z, K, wc, wb):
	'''
	transfer function:
	    1+(s/wb)	1
	 Z=K-------- --------
	     (s/wb)	 1+(s/wc)				'''
	Z_model = K*tf.lti([1/wb,1], [1/wb,0])*tf.lti([1], [1/wc,1])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def Serial_RC_freq(W_z, K, wc):
	'''Model corresponding to:
	--R--C--
	with a transfer function:
	    1+(s/wb)
	 Z=K--------
	     (s/wb)					'''
	Z_model = K*tf.lti([1/wc,1], [1/wc,0])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def Serial_RC_freq2(W_z, K, wc):
	'''Model corresponding to:
	--R--C--
	with a transfer function:
	    1+(s/wb)
	 Z=K--------
	     (s/wb)					'''
	Z_model = K*tf.lti([1/wc,1], [1/wc,0])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def Serial_RCPE_freq(W_z, gamma, K, wb):
	'''Model corresponding to:
	--R--CPE--
	with a transfer function:
	    1+(s/wb)**gamma
	 Z=K---------------
	     (s/wb)**gamma			'''
	D = ne.oustaloup_approx(gamma,Wl,Wh,N_zp)
	Z_model = K*(((D/(wb**gamma))+1)/(D/(wb**gamma)))
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def Serial_RCPE_freq2(W_z, gamma, K, wb):
	'''Model corresponding to:
	--R--CPE--
	with a transfer function:
	    1+(s/wb)**gamma
	 Z=K---------------
	     (s/wb)**gamma			'''
	D = ne.oustaloup_approx(gamma,Wl,Wh,N_zp)
	Z_model = K*(((D/(wb**gamma))+1)/(D/(wb**gamma)))
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def RCR_Non_fractional_freq(W_z,K,wn,wd):
	'''Model corresponding to:
	     |--R--|
	--R--|	   |--
	     |--C--|
	with a transfer function:
	    1+(s/wn)
	 Z=K---------
	    1+(s/wd)				'''
	Z_model = K*tf.lti([1/wn,1], [1/wd,1])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def RCR_Non_fractional_freq2(W_z,K,wn,wd):
	'''Model corresponding to:
	     |--R--|
	--R--|	   |--
	     |--C--|
	with a transfer function:
	    1+(s/wn)
	 Z=K---------
	    1+(s/wd)				'''
	Z_model = K*tf.lti([1/wn,1], [1/wd,1])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def Serial_QQR_freq(W_z, gamma_1, wb_1, gamma_2, wb_2, R):
	'''Model corresponding to:
	--R--CPE--CPE---
	'''
	D_1 = ne.oustaloup_approx(gamma_1,Wl,Wh,N_zp)
	D_2 = ne.oustaloup_approx(gamma_2,Wl,Wh,N_zp)
	Z_model = 1/(D_1/(wb_1**gamma_1))+1/(D_2/(wb_2**gamma_2))+R
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def Serial_QQR_freq2(W_z, gamma_1, wb_1, gamma_2, wb_2, R):
	'''Model corresponding to:
	--R--CPE--CPE---
	'''
	D_1 = ne.oustaloup_approx(gamma_1,Wl,Wh,N_zp)
	D_2 = ne.oustaloup_approx(gamma_2,Wl,Wh,N_zp)
	Z_model = 1/(D_1/(wb_1**gamma_1))+1/(D_2/(wb_2**gamma_2))+R
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def Low_pass_freq(W_z, K, wc):
	'''Model corresponding to:
	--R--C--
	with a transfer function:
	       1
	 Z=K--------
	    1+(s/wb)					'''
	Z_model = K*tf.lti([1], [1/wc,1])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def Low_pass_freq2(W_z, K, wc):
	'''Model corresponding to:
	--R--C--
	with a transfer function:
	       1
	 Z=K--------
	    1+(s/wb)					'''
	Z_model = K*tf.lti([1], [1/wc,1])
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def Fractional_low_pass_freq(W_z, K, wc, gamma):
	'''Model corresponding to:
	--R--C--
	with a transfer function:
	       1
	 Z=K--------
	    1+(s/wb)					'''
	D = ne.oustaloup_approx(gamma,Wl,Wh,N_zp)
	#Z_model = K*tf.lti([1], [1/wc,1])
	Z_model = K*((1.)/(1+(D/(wc**gamma))))
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return mag_z

def Fractional_low_pass_freq2(W_z, K, wc, gamma):
	'''Model corresponding to:
	--R--C--
	with a transfer function:
	       1
	 Z=K--------
	    1+(s/wb)					'''
	D = ne.oustaloup_approx(gamma,Wl,Wh,N_zp)
	#Z_model = K*tf.lti([1], [1/wc,1])
	Z_model = K*((1.)/(1+(D/(wc**gamma))))
	w_z, mag_z, phase_z = tf.bode(Z_model,w=W_z)
	return phase_z

def RCPE_pole_freq(W_z, gamma, K, wb, wc):
	''' RCPE model with a lowpass behaviour in HF:
	       |--R--|
	--CPE--|	 |--
	       |--C--|			'''
	D = ne.oustaloup_approx(gamma,Wl,Wh,N_zp)
	Z_1 = K*(((D/(wb**gamma))+1)/(D/(wb**gamma)))
	Z_2 = tf.lti([1], [1/wc,1])
	w_z, mag_1, phase_1 = tf.bode(Z_1,w=W_z)
	w_z, mag_2, phase_2 = tf.bode(Z_1,w=W_z)
	mag_z = mag_1 + mag_2
	return mag_z

def RCPE_pole_freq2(W_z, gamma, K, wb, wc):
	''' RCPE model with a lowpass behaviour in HF:
	       |--R--|
	--CPE--|	 |--
	       |--C--|			'''
	D = ne.oustaloup_approx(gamma,Wl,Wh,N_zp)
	Z_1 = K*(((D/(wb**gamma))+1)/(D/(wb**gamma)))
	Z_2 = tf.lti([1], [1/wc,1])
	w_z, mag_1, phase_1 = tf.bode(Z_1,w=W_z)
	w_z, mag_2, phase_2 = tf.bode(Z_1,w=W_z)
	phase_z = phase_1 + phase_2
	return phase_z

def fractional_alpha_beta_freq(W_z, gamma1, K, wb, wc, gamma2 ):
	''' Model with fractional low frequency 
	and low pass filter behaviour - 
	cf work with A. Degache - IMS'''
	mod_1 = Serial_RCPE_freq(W_z, gamma1, K, wb)
	mod_2 = Fractional_low_pass_freq(W_z, 1.0, wc, gamma2)
	return mod_1+mod_2

def fractional_alpha_beta_freq2(W_z, gamma1, K, wb, wc, gamma2 ):
	''' Model with fractional low frequency 
	and low pass filter behaviour - 
	cf work with A. Degache - IMS'''
	phase_1 = Serial_RCPE_freq2(W_z, gamma1, K, wb)
	phase_2 = Fractional_low_pass_freq2(W_z, 1.0, wc, gamma2)
	return phase_1+phase_2
