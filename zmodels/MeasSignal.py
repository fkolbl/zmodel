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
# Last changed: 2018-06-13														#
#################################################################################
import numpy as np
import scipy.optimize as optim
import zmodels.Models as mod
import zmodels.norm as norm
import logging

# Constants
epsilon_0 = 4*np.pi*10e-7

## LOGGING FILE HANDLING
verbosity = 'high'
loggname = 'zmodel_activity.log'

if verbosity == 'high':
	logging.basicConfig(filename=loggname,level=logging.DEBUG)
if verbosity=='medium':
	logging.basicConfig(filename=loggname,level=logging.INFO)
else:
	logging.basicConfig(filename=loggname)

## CLASS
class MeasSignal(object):
	"""docstring for MeasSignal"""
	def __init__(self, filename, meas_type,comment=None):
		super(MeasSignal, self).__init__()
		self.source = filename
		self.meas_type = meas_type

		# testing that the filename corresponds to an actual file and opening it
		if not('.csv' in self.source):
			logging.warning('Warning: Measurements results should be in comma separated values files')
		try:
			test = open(self.source)
		except IOError:
			logging.error('Error:  file %s not found'%self.source)
			quit()
		else:
			del test
			self.data = np.genfromtxt(self.source, delimiter=',')
			# Fields for measurement data collection
			self.freq = None
			self.zimp = None
			self.zimp_dB = None
			self.zphase = None
			self.zRe = None
			self.zIm = None
			self.yphase = None
			self.yimp = None
			self.yRe = None
			self.yIm = None
			self.sigma = None
			self.epsilon_r = None
			self.I_stim = None
			self.V_stim = None
			self.I_meas = None
			self.V_meas = None
			self.time = None
			self.comment = comment
			self.electrode_Area = None
			self.electrode_Distance = None
			if self.meas_type=='freq':
				# Frequential measurement
				self.freq = self.data[:,0]
				self.zimp = self.data[:,1]
				self.zimp_dB = 20*np.log10(self.zimp)
				if len(self.data[0])>2:
					# computation of the real/imaginary part of the impedance
					self.zphase = self.data[:,2]
					self.zRe = self.zimp*np.cos(self.zphase*(np.pi/180.))
					self.zIm = self.zimp*np.sin(self.zphase*(np.pi/180.))
					# computation of the complex admittance
					self.yphase = -self.zphase
					self.yimp = 1/self.zimp
					self.yRe = self.yimp*np.cos(self.yphase*(np.pi/180.))
					self.yIm = self.yimp*np.sin(self.yphase*(np.pi/180.))
				else:
					self.zphase = []
			elif sel.meas_type=='temp':
				# Temporal measurement with voltage solicitation
				logging.debug('not implemented yet')
			else:
				logging.error('Error: unrecognized measurement type, operation aborted')
				quit()
		# cleaning memory
		del self.data
		# fields for the identification/parameter estimation process
		self.ident = False
		self.model = []
		self.method = []
		self.param = []
		self.cov = []
		self.zimp_dB_ident = []
		self.zimp_ident = []
		self.zphase_ident = []
		self.V_ident = []
		self.I_ident = []
		self.nrmse_models = []

	def compute_permittivity(self):
		if (self.electrode_Distance == None or self.electrode_Area==None):
			logging.error('Error: Electrode geometry definition (Area and Distance) requiered for permittity computation')
			quit()
		else:
			self.geom_coeff = self.electrode_Area/self.electrode_Distance
			self.sigma = self.yRe/self.geom_coeff
			self.epsilon_r = (self.yRe/(self.geom_coeff*2*np.pi*self.freq*epsilon_0))

	def identif(self,modelname,method,p0=None,bounds=(-np.inf, np.inf),sub_method=None):
		if method == 'LeastSquare':
			return self.least_square_fit(modelname,p0=p0,bounds=bounds,method=sub_method)
		elif method == 'MinimizeError':
			if bounds==(-np.inf, np.inf):
				return self.minimize_error_fit(modelname,p0 = p0,method=sub_method)
			else:
				print bounds
				return self.minimize_error_fit(modelname,p0 = p0, bounds = bounds,method=sub_method)
		else:
			logging.error('Error: unsupported identification/parameter estimation method')
			quit()

	def least_square_fit(self,modelname,p0=None,bounds=(-np.inf, np.inf),method=None,quantity='mod'):
		if self.meas_type=='freq':
			# load the correct model function
			func_name = modelname+'_'+self.meas_type
			used_model = getattr(mod, func_name)
			model_phase = getattr(mod, modelname+'_freq2')
			popt, pcov = optim.curve_fit(used_model, self.freq, self.zimp_dB, p0=p0, bounds=bounds, method=method)#
			# store results
			self.ident = True
			self.model.append(modelname)
			self.method.append('LeastSquare')
			self.param.append(popt)
			self.cov.append(pcov)
			self.zimp_dB_ident.append(used_model(self.freq,*popt))
			self.zimp_ident.append(10**(used_model(self.freq,*popt)/20))
			self.zphase_ident.append(model_phase(self.freq,*popt))
			# saved results in logg
			logging.info('Least Square Fit on file %s'%(self.source))
			logging.info('\t found for model %s the parameters: %s'%(self.model[-1],str(self.param[-1])))
			logging.info('\t\t the Max Error is : %s'%norm.Max_Error(self.zimp,self.zimp_ident[-1]))
			logging.info('\t\t the Mean Error is : %s'%norm.Mean_Error(self.zimp,self.zimp_ident[-1]))
			logging.info('\t\t the MSE is : %s'%norm.MSE(self.zimp,self.zimp_ident[-1]))
			logging.info('\t\t the RMSE is : %s'%norm.RMSE(self.zimp,self.zimp_ident[-1]))
			logging.info('\t\t the NRMSE is : %s'%norm.NRMSE(self.zimp,self.zimp_ident[-1]))

			# return parameters and covariance
			return popt, pcov
		elif self.meas_type=='temp_volt':
			logging.debug('not implemented yet')
			quit()
		else:
			logging.debug('not implemented yet')
			quit()
		return pcov

	def minimize_error_fit(self,modelname,p0=None,bounds=None,method=None,quantity='mod'):
		if self.meas_type=='freq':
			if quantity=='phase':
				func_name = modelname+'_'+self.meas_type+'2'
			else:
				func_name = modelname+'_'+self.meas_type
			print func_name
			used_model = getattr(mod, func_name)
			model_phase = getattr(mod, modelname+'_freq2')
			fmin = self.define_error_function(modelname)
			res = optim.minimize(fmin, list(p0), options={'disp': True},bounds=bounds)
			# store results
			self.ident = True
			self.model.append(modelname)
			self.method.append('MinimizeError')
			self.param.append(res.x)
			self.cov.append(None)
			self.zimp_dB_ident.append(used_model(self.freq,*res.x))
			self.zimp_ident.append(10**(used_model(self.freq,*res.x)/20))
			self.zphase_ident.append(model_phase(self.freq,*res.x))
			return res.x
		elif self.meas_type=='temp_volt':
			logging.debug('not implemented yet')
			quit()
		else:
			logging.debug('not implemented yet')
			quit()

	def define_error_function(self,modelname):
		func_name = modelname+'_'+self.meas_type
		used_model = getattr(mod, func_name)
		def fmin(x):
			mag_model = used_model(self.freq,*x)
			nrmse = norm.NRMSE(self.zimp_dB,mag_model)
			return nrmse
		return fmin
