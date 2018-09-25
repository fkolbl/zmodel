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
from scipy import signal
from numpy import polymul,polyadd

class lti(object):
	"""Linear time invariant class,
	This class is directly based on the SCIPY LTI class. 
	Unless the SCIPY version, algebraic operation are enabled, however most of the functions
	have the same name and parameters than their SCIPY implementation.
	Temporal and frequential resolutions are also performed using directly SCIPY"""
	def __init__(self, *arg):
		super(lti, self).__init__()
		self.arg = arg
		# preventing errors
		self.num = 0
		self.den = 0
		# computation of the LTI system
		if len(arg) == 1:
			print 'error: a lti needs at least 2 arguments to be defined'
		elif len(arg) > 3:
			print 'error: too many arguments'
		else:
			# No error, if two arguments, defining as numerator, denominator
			if len(arg) == 2:
				self.scipylti = signal.lti(self.arg[0], self.arg[1])
				self.num = self.scipylti.num
				self.den = self.scipylti.den
			# if three arguments, defining as zero, pole gain
			if len(arg) == 3:
				self.scipylti = signal.lti(self.arg[0], self.arg[1], self.arg[2])
				self.num = self.scipylti.num
				self.den = self.scipylti.den
		# cleaning the memory
		del self.arg
		del self.scipylti

	# Algebric operations on LTI redefinintion
	def __neg__(self):
		return lti(-self.num,self.den)

	def __add__(self,other):
		if type(other) in [type(self)]:
			numer = polyadd(polymul(self.num,other.den),polymul(other.den,self.num))
			denom = polymul(self.den,other.den)
			return lti(numer,denom)
		else:
			return lti(polyadd(self.num,self.den*float(other)),self.den)
		

	__radd__ = __add__ # symmetric behaviour as the addition is commutative

	def __sub__(self,other):
		if type(other) in [type(self)]:
			numer = polyadd(polymul(self.num,other.den),-polymul(other.den,self.num))
			denom = polymul(self.den,other.den)
			return lti(numer,denom)
		else:
			return lti(polyadd(self.num,-self.den*float(other)),self.den)

	def __rsub__(self,other):
		if type(other) in [type(self)]:
			numer = polyadd(polymul(other.num,self.den),-polymul(self.den,other.num))
			denom = polymul(self.den,other.den)
			return lti(numer,denom)
		else:
			return lti(polyadd(-self.num,self.den*float(other)),self.den)

	def __mul__(self,other):
		if type(other) in [type(self)]:
			numer = polymul(self.num,other.num)
			denom = polymul(self.den,other.den)
			return lti(numer,denom)
		if type(other) in [int, float]:
			return lti(self.num*float(other),self.den)

	__rmul__ = __mul__ # symmetric behaviour as the multiplication is commutative

	def __div__(self,other):
		if type(other) in [type(self)]:
			numer = polymul(self.num,other.den)
			denom = polymul(self.den,other.num)
			return lti(numer,denom)
		else:
			return lti(self.num,self.den*float(other))

	def __rdiv__(self,other):
		if type(other) in [type(self)]:
			numer = polymul(self.den,other.num)
			denom = polymul(self.num,other.den)
			return ltimul(numer,denom)
		else:
			return lti(float(other)*self.den,self.num)

	# Method to recompute Scipy LTI form to re-hinerit operations and methods from Scipy
	def pass2scipy(self):
		self.hidden_lti = signal.lti(self.num, self.den)

	# Invert Method of the pass2scipy
	def scipy2lti(self):
		self.num = self.hidden_lti.num
		self.den = self.hidden_lti.den

	# Clean up the mess
	def clean_hidden(self):
		del self.hidden_lti

# Functions for LTIs
def bode(LTI, w=None, n=100):
	'''Clone function of the scipy.signal bode function,
	adapted to LTI with enabled algebraic operations,
	returns the bode of a LTI, parameters are:
		LTI: an instance of the zmodels LTI class
		w: array_like, optional,
			Array of frequencies (in rad/s).
			Magnitude and phase data is calculated for every value in this array.
			If not given a reasonable set will be calculated.
		n: int, optional
			Number of frequency points to compute if w is not given.
			The n frequencies are logarithmically spaced in an interval chosen to 
			include the influence of the poles and zeros of the system.
	outputs are:
		w: 1D ndarray
			Frequency array [rad/s]
		mag: 1D ndarray
			Magnitude array [dB] 
		phase: 1D ndarray
			Phase array [deg]
	'''
	LTI.pass2scipy()
	w, mag, phase = signal.bode(LTI.hidden_lti, w=w, n=n)
	LTI.clean_hidden
	return w, mag, phase

def nyquist(LTI, w=None, n=100):
	'''Returns the real and imaginary parts of a LTI frequencial response,
	parameters are:
		LTI: an instance of the zmodels LTI class
		w: array_like, optional,
			Array of frequencies (in rad/s).
			Magnitude and phase data is calculated for every value in this array.
			If not given a reasonable set will be calculated.
		n: int, optional
			Number of frequency points to compute if w is not given.
			The n frequencies are logarithmically spaced in an interval chosen to 
			include the influence of the poles and zeros of the system.
	outputs are:
		w: 1D ndarray
			Frequency array [rad/s]
		Re: 1D ndarray
			Real part array 
		phase: 1D ndarray
			Imaginary part array
	'''
	LTI.pass2scipy()
	w, mag, phase = signal.bode(LTI.hidden_lti, w=w, n=n)
	LTI.clean_hidden
	mod = Magnitude2Modulus(mag)
	Re, Im = Polar2Algebraic(mod, phase)
	return w, Re, Im

def impulse(LTI, X0=None, T=None, N=None):
	'''Clone function of the scipy.signal impulse function,
	adapted to LTI with enabled algebraic operations,
	returns the impulse response of a LTI. parameters are:
		LTI: an instance of the zmoel LTI class
		X0: array_like, optional
			Initial state-vector (default is zero)
		T: array_like, optional
			Time points (computed if not given).
		N: int, optional
			Number of time points to compute if T is not given.
	Outputs are:
		t: 1D ndarray
			Output time points.
		yout: 1D ndarray
			Impulse response of system.
	'''
	LTI.pass2scipy()
	t, yout = signal.impulse(LTI.hidden_lti, X0=X0, T=T, N=N)
	LTI.clean_hidden
	return t, yout

def impulse2(LTI, X0=None, T=None, N=None, **kwargs):
	'''Clone function of the scipy.signal impulse2 function,
	adapted to LTI with enabled algebraic operations,
	returns the impulse response of a LTI. parameters are:
		LTI: an instance of the zmoel LTI class
		X0: array_like, optional
			Initial state-vector (default is zero)
		T: array_like, optional
			Time points (computed if not given).
		N: int, optional
			Number of time points to compute if T is not given.
		kwargs : various types
			Additional keyword arguments are passed on the function scipy.signal.lsim2,
 			see scipy documentation for more details 
	Outputs are:
		t: 1D ndarray
			Output time points.
		yout: 1D ndarray
			Impulse response of system.
	'''
	LTI.pass2scipy()
	t, yout = signal.impulse2(LTI.hidden_lti, X0=X0, T=T, N=N,**kwargs)
	LTI.clean_hidden
	return t, yout

def lsim(LTI, U, T, X0=None, interp=True):
	LTI.pass2scipy()
	T, yout, xout = signal.lsim(LTI.hidden_lti, U, T, X0=X0, interp=interp)
	LTI.clean_hidden
	return T, yout, xout

def lsim2(LTI, U=None, T=None, X0=None, **kwargs):
	LTI.pass2scipy()
	T, yout, xout = signal.lsim2(LTI.hidden_lti, U=U, T=T, X0=X0,**kwargs)
	LTI.clean_hidden
	return T, yout, xout	

def step(LTI, X0=None, T=None, N=None):
	'''Clone function of the scipy.signal step function,
	adapted to LTI with enabled algebraic operations,
	returns the step response of a LTI. parameters are:
		LTI: an instance of the zmoel LTI class
		X0: array_like, optional
			Initial state-vector (default is zero)
		T: array_like, optional
			Time points (computed if not given).
		N: int, optional
			Number of time points to compute if T is not given.
	Outputs are:
		t: 1D ndarray
			Output time points.
		yout: 1D ndarray
			Step response of system.
	'''
	LTI.pass2scipy()
	t, yout = signal.step(LTI.hidden_lti, X0=X0, T=T, N=N)
	LTI.clean_hidden
	return t, yout

def step2(LTI, X0=None, T=None, N=None, **kwargs):
	'''Clone function of the scipy.signal step2 function,
	adapted to LTI with enabled algebraic operations,
	returns the step response of a LTI. parameters are:
		LTI: an instance of the zmoel LTI class
		X0: array_like, optional
			Initial state-vector (default is zero)
		T: array_like, optional
			Time points (computed if not given).
		N: int, optional
			Number of time points to compute if T is not given.
		kwargs : various types
			Additional keyword arguments are passed on the function scipy.signal.lsim2,
 			see scipy documentation for more details 
	Outputs are:
		t: 1D ndarray
			Output time points.
		yout: 1D ndarray
			Step response of system.
	'''
	LTI.pass2scipy()
	t, yout = signal.step2(LTI.hidden_lti, X0=X0, T=T, N=N, **kwargs)
	LTI.clean_hidden
	return t, yout

def Magnitude2Modulus(Magnitude):
	if type(Magnitude) in [float,int]:
		mod = 10**(Magnitude/20)
	else:
		mod = np.zeros(len(Magnitude))
		for  k in xrange(len(Magnitude)):
			mod[k] = 10**(Magnitude[k]/20)
	return mod

def Polar2Algebraic(Modulus, Phase):
	if type(Modulus) in [float,int]:
		Re = Modulus*math.cos(Phase*math.pi/180)
		Im = Modulus*math.sin(Phase*math.pi/180)
	else:
		Re = np.zeros(len(Modulus))
		Im = np.zeros(len(Modulus))
		for  k in xrange(len(Modulus)):
			Re[k] = Modulus[k]*math.cos(Phase[k]*math.pi/180)
			Im[k] = Modulus[k]*math.sin(Phase[k]*math.pi/180)
	return Re, Im
