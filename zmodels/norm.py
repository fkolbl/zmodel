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
# Last changed: 2017-03-01														#
#################################################################################
import math
import numpy as np

def N_inf(vector):
	return np.amax(vector)

def Norm(vector, ind):
	if ind == 'inf':
		return N_inf(vector)
	elif ind<=0:
		print 'norms defined for strictely positive indices'
		quit()
	else:
		temp = 0
		for x in vector:
			temp += abs(x)**ind
		return temp**(1/float(ind))

def Eucl_Norm(vector):
	return Norm(vector,2.0)

def dist(P1, P2):
	if len(P1)!=len(P2):
		print 'distance computation error: dimension mismatch'
		quit()
	else:
		diff = P1 - P2
		return Eucl_Norm(diff)

def Max_Error(signal1,signal2):
	error = abs(signal2 - signal1)
	return np.max(error)

def Mean_Error(signal1,signal2):
	error = signal2 - signal1
	return abs(np.sum(error)/len(error))

def MSE(signal1,signal2):
	error = signal2 - signal1
	return np.sum(np.square(error))/len(error)


def RMSE(signal1, signal2):
	error = signal2 - signal1
	return np.sqrt(MSE(signal1,signal2))

def NRMSE(signal1, signal2):
	error = signal2 - signal1
	normalisation = np.max(signal1) - np.min(signal1)
	return np.sqrt(MSE(signal1,signal2))/normalisation

