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
# Last changed: 2016-11-04														#
#################################################################################
import math
import numpy as np
import zmodels.lti as tf
import sympy as sp
z = sp.symbols('z')
sp.init_printing(use_unicode=False, wrap_line=False, no_global=True)

def theoretical_bode(gamma, w, Wu = 1):
	# Returns the magnitude and phase of a ideal fractional derivator as:
	# mag, phase
	# Arguments are:
	#	gamma	: the fractional integration order
	#	w 		: a frequency vector on which computation are conduced
	#	Wu		: optional, corresponds to the unity gain frequency
	#			  if not specified, corresponds to 1 rad/s

	mag = np.zeros(len(w))
	phase = np.zeros(len(w))
	for k in range(len(w)):
		# print (w[k]/Wu)**gamma
		mag[k] = 20*math.log10((w[k]/Wu)**gamma)
		phase[k] = 90*gamma
	return mag, phase

def oustaloup_approx(gamma,Wb,Wh,Nzp=0,unit = True):
	# Returns the Ousltaloup approximation of a fractional derivator
	# as a Scipy LTI system.
	# Arguments are:
	# 	gamma	: the fractional derivation order
	#	Wb		: the lowest frequency of the approximation domain (in rad/s)
	#	Wh		: the highest frequency of the approximation domain (in rad/s)
	# 	Nzp		: the number of zeros/poles on the approximation, 
	# 			  note that this parameter is optional, if not specified there are 
	#			  two pairs ot zeros/poles per decade
	# 	unit	: if unit = True, the unity gain is for omega = 0
	# 				 unit = False, the unity fain is for the middle frequency band between
	#				     Wb and Wh

	# automatic assignement of 2 zeros/poles per decade if not specified
	if Nzp == 0:
		Nzp = int(2*math.ceil(math.log10(Wh) - math.log10(Wb)))
	# quantities used in the Oustaloup approximation
	alpha_eta=(Wh/Wb)**(1.0/Nzp);
	alpha=alpha_eta**gamma;
	eta=alpha_eta/alpha;
	delta = np.zeros(Nzp)
	for k in range(len(delta)):
		delta[k] = math.sqrt(eta)*alpha_eta**(k)
	# Computing the zeros and poles
	Wzero = delta*Wb
	Wpole = Wh/delta
	# Computing the static gain, discriminating unity gain at the middle frequency or at wu=1rad/s
	Co = (Wb/math.sqrt(Wb*Wh))**gamma
	if unit:
		Co = (Co/math.sqrt(Wb*Wh)**gamma)
	return tf.lti(-Wzero, -Wpole, 1/Co)

def A(z,gamma,n):
	if n == 0:
		A0 = 1
		return A0
	else:
		if n%2 == 0:
			cn = 0
		else:
			cn = float(gamma)/n
		An = A(z,gamma,n-1)-cn*z**(-n)*A(z**(-1),gamma,n-1)
		return sp.expand(An)

def recursiv_expanded_Tustin_Poly(gamma,n):
	return sp.expand(A(z,gamma,n)*(z**n))

def recursiv_expanded_Tustin_Coeff(gamma,n):
	P = sp.Poly(recursiv_expanded_Tustin_Poly(gamma,n))
	coef_list = P.coeffs()
	newP = []
	for item in coef_list:
		newP.append(float(item))
	return np.array(newP)

def expandedTustin_approx(gamma,N,Tech):
	b = ((2./Tech)**gamma)*recursiv_expanded_Tustin_Coeff(gamma,N)
	a = recursiv_expanded_Tustin_Coeff(-gamma,N)
	return b,a