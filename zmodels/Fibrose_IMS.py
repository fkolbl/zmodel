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
# First added:  2017-12-12														#
# Last changed: 2018-06-13														#
#################################################################################
import zmodels.MeasSignal as expe
import zmodels.Models as mod
import zmodels.lti as tf
import zmodels.fractional as ne
import zmodels.norm as nrm

import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import sys
import os
import glob
import math
import pandas as pd

### FUNCTIONS TO HANDLE FILE FORMAT AND CONVENTIONS USED AT IMS
def z2csv(zfilename,csvfilename):
	with open(zfilename,'r') as infile, open(csvfilename, 'w') as outfile:
		for line in infile:
			word_list = line.split()
			try:
				float(word_list[0])
				freq_str = word_list[0]
				Re_str = word_list[4]
				Im_str = word_list[5]
				Re = float(Re_str)
				Im = float(Im_str)
				Mod = math.hypot(Re,Im)
				Phi = math.degrees(math.atan2(Im,Re))
				line_to_write = freq_str+','+str(Mod)+','+str(Phi)+'\n'
				outfile.write(line_to_write)
			except ValueError:
				pass

def check_file_name(filename,path=''):
	list_of_items = filename.replace(path,'').replace('/','').replace('.csv','').split('_')
	date = list_of_items[0]
	electrode = list_of_items[1]
	amplitude = 0
	if 'uV' in list_of_items[2]:
		amplitude = float(list_of_items[2].replace('uV',''))*1e-6
	elif 'mV' in list_of_items[2]:
		amplitude = float(list_of_items[2].replace('mV',''))*1e-3
	elif 'V' in list_of_items[2]:
		amplitude = float(list_of_items[2].replace('V',''))
	tissu = list_of_items[4]
	return [date,electrode,amplitude,tissu]

### CLASS 
class ExVivo_Experiment(object):
	"""docstring for ExVivo_Experiment"""
	def __init__(self, Path, Folders, out_name='zmodels_out'):
		super(ExVivo_Experiment, self).__init__()
		self.Path = Path
		self.Folders = Folders

		self.File_list = []
		self.Data_list = []

		# files to save results
		self.outfile_RC_N = out_name+'_Serial_RC.csv'
		self.outfile_RCPE_N = out_name+'_Serial_RCPE.csv'
		self.outfile_RCZeroPole_N = out_name+'_RC_Zero_Pole.csv'
		self.outfile_fractional_alpha_beta_N = out_name+'_fractional_alpha_beta.csv'
		self.outfile_bestfit_N = out_name+'_bestfit.csv'
		# create or empty file, write headers and close for later use
		self.outfile_RC_F = open(self.outfile_RC_N,'w')
		self.outfile_RCPE_F = open(self.outfile_RCPE_N,'w')
		self.outfile_RCZeroPole_F = open(self.outfile_RCZeroPole_N,'w')
		self.outfile_fractional_alpha_beta_F = open(self.outfile_fractional_alpha_beta_N,'w')
		self.outfile_bestfit_F = open(self.outfile_bestfit_N,'w')
		self.outfile_RC_F.write('source\t K\t wc\t NRMSE \n')
		self.outfile_RCPE_F.write('source\t gamma\tK\t wc\t NRMSE  \n')
		self.outfile_RCZeroPole_F.write('source\t K\t wc\t wb\t NRMSE \n')
		self.outfile_fractional_alpha_beta_F.write('source\t gamma1\t K\t wb\t wc\t gamma2\t NRMSE\n')
		self.outfile_bestfit_F.write('source\t model\t NRMSE\n')
		self.outfile_RC_F.close()
		self.outfile_RCPE_F.close()
		self.outfile_RCZeroPole_F.close()
		self.outfile_fractional_alpha_beta_F.close()
		self.outfile_bestfit_F.close()
		# meta data
		self.Date_list = []
		self.Electrode_list = []
		self.Amplitude_list = []
		self.Tissu_list = []

	def load_all_measurements(self):
		# convert zfiles to csv to adapt to MeasSignal class and load the corresponding data
		print 'Converting all z files to csv files'
		for folder in self.Folders:
			print '\t'+self.Path+folder
			for infile in glob.glob( os.path.join(self.Path+folder, '*.z') ):
				csv_file_path = infile.replace('.z','.csv')
				z2csv(infile,csv_file_path)
				self.File_list.append(csv_file_path)
				print "loading %s"%csv_file_path	
				self.Data_list.append(expe.MeasSignal(csv_file_path,'freq'))
				meta_data = check_file_name(csv_file_path,path = self.Path+folder)
				self.Date_list.append(meta_data[0])
				self.Electrode_list.append(meta_data[1])
				self.Amplitude_list.append(meta_data[2])
		self.Data_NRMSE = np.zeros(len(self.Data_list))

	def plot_all_measurements(self,pict_format ='.pdf'):
		for data in self.Data_list:
			print "saving %s"%data.source.replace('.csv',pict_format)
			plt.figure()
			plt.loglog(data.freq, data.zimp)
			plt.grid(True, which="both")
			plt.xlabel('frequency (Hz)')
			plt.ylabel('Impedance magnitude')
			plt.savefig(data.source.replace('.csv',pict_format))
			plt.close()

	def identify_models(self):
		self.outfile_RC_F = open(self.outfile_RC_N,'a')
		self.outfile_RCPE_F = open(self.outfile_RCPE_N,'a')
		self.outfile_RCZeroPole_F = open(self.outfile_RCZeroPole_N,'a')
		self.outfile_fractional_alpha_beta_F = open(self.outfile_fractional_alpha_beta_N,'a')
		self.outfile_bestfit_F = open(self.outfile_bestfit_N,'a')
		for data in self.Data_list:
			print "... identifying on file %s"%data.source
			try:
				# Simple one zero model for initial conditions
				data.least_square_fit('Serial_RC',bounds=([0.0, 0.0],[np.inf, np.inf]))
				data.nrmse_models.append(nrm.NRMSE(data.zimp,data.zimp_ident[-1]))
				# dynamic initial guess for RCPE model
				gamma_dyn = 1
				K_dyn = data.param[-1][0]
				wb_dyn = data.param[-1][1]
				self.outfile_RC_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\n')
				try:
					data.least_square_fit('Serial_RCPE', p0=(gamma_dyn, K_dyn, wb_dyn))
					data.nrmse_models.append(nrm.NRMSE(data.zimp,data.zimp_ident[-1]))
					self.outfile_RCPE_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\n')
				except RuntimeError:
					print 'identification not successful using Least Square algorithm for simple RC model'
					print '... relaunching with gradient descend'
					try:
						data.minimize_error_fit('Serial_RCPE', p0=(gamma_dyn, K_dyn, wb_dyn))
						data.nrmse_models.append(nrm.NRMSE(data.zimp,data.zimp_ident[-1]))
						self.outfile_RCPE_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\n')
					except RuntimeError:
						print '... failed to identify Serial_RCPE model'
						self.outfile_RCPE_F.write(data.source+'\t 0.0 \t 0.0 \t 0.0 \n')
			except RuntimeError:
				print 'identification not successful using Least Square algorithm for simple RC model'
				self.outfile_RC_F.write(data.source+'\t 0.0 \t 0.0 \n')
				pass
			try:
				# Simple one zero one pode for initial conditions
				data.least_square_fit('RC_Zero_Pole',bounds=([0.0, 0.0, 0.0],[np.inf, np.inf, np.inf]))
				data.nrmse_models.append(nrm.NRMSE(data.zimp,data.zimp_ident[-1]))
				self.outfile_RCZeroPole_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\n')
				# dynamic initial guess for complex model
				gamma1_dyn = 0.9
				K_dyn = data.param[-1][0]
				wb_dyn = data.param[-1][2]
				wc_dyn = data.param[-1][1]
				gamma2_dyn = 0.5
				try:
					data.least_square_fit('fractional_alpha_beta', p0=(gamma1_dyn, K_dyn, wb_dyn, wc_dyn,gamma2_dyn))
					data.nrmse_models.append(nrm.NRMSE(data.zimp,data.zimp_ident[-1]))
					self.outfile_fractional_alpha_beta_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\t'+str(data.param[-1][3])+'\t'+str(data.param[-1][4])+'\n')
				except RuntimeError:
					print 'indentification not successful using Least Square algorithm'
					print '... relaunching with gradient descend'
					try:
						data.minimize_error_fit('fractional_alpha_beta', p0=(gamma1_dyn, K_dyn, wb_dyn, wc_dyn,gamma2_dyn))
						data.nrmse_models.append(nrm.NRMSE(data.zimp,data.zimp_ident[-1]))
						self.outfile_fractional_alpha_beta_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\t'+str(data.param[-1][3])+'\t'+str(data.param[-1][4])+'\n')
					except RuntimeError:
						print '... failed to identify fractional_alpha_beta model'
						self.outfile_fractional_alpha_beta_F.write(data.source+'\t 0.0 \t 0.0\t 0.0\t 0.0\n')
			except RuntimeError:
				self.outfile_RCZeroPole_F.write(data.source+'\t 0.0 \t 0.0\t 0.0\n')
				print 'identification not successful using Least Square algorithm for RC Zero and Pole model'

			# check for the model with the best fit
			models_nrmse = np.zeros(len(data.zimp_ident))
			for k in xrange(len(data.zimp_ident)):
				models_nrmse[k] = nrm.NRMSE(data.zimp,data.zimp_ident[k])
			minNRMSEarg = np.argmin(models_nrmse)
			print models_nrmse
			print np.asarray(data.nrmse_models)
			print str(models_nrmse[minNRMSEarg])
			self.outfile_bestfit_F.write(data.source+'\t' + data.model[minNRMSEarg]+'\t' + str(models_nrmse[minNRMSEarg]) +'\n')
			self.Data_NRMSE[self.Data_list.index(data)] = models_nrmse[minNRMSEarg]
			
			# plot results
			plt.figure()
			plt.loglog(data.freq, data.zimp, label='measurement')
			for k in xrange(len(data.zimp_ident)):
				plt.loglog(data.freq, data.zimp_ident[k], label=data.model[k])
			plt.grid(True, which="both")
			plt.legend()
			plt.xlabel('frequency (Hz)')
			plt.ylabel('Impedance magnitude')
			plt.savefig(data.source.replace('.csv','_indent_mod.pdf'))
			plt.close()

			plt.figure()
			plt.semilogx(data.freq, data.zphase, label='measurement')
			for k in xrange(len(data.zphase_ident)):
				plt.semilogx(data.freq, data.zphase_ident[k], label=data.model[k])
			plt.grid(True, which="both")
			plt.legend(loc=0)
			plt.xlabel('frequency (Hz)')
			plt.ylabel('Impedance phase')
			plt.savefig(data.source.replace('.csv','_indent_phase.pdf'))
			plt.close()
		self.outfile_RC_F.close()
		self.outfile_RCPE_F.close()
		self.outfile_RCZeroPole_F.close()
		self.outfile_fractional_alpha_beta_F.close()
		self.outfile_bestfit_F.close()
		print 'NRMSE results for the set of data used: '
		print 'The maximum NRMSE value found was: '+str(np.max(self.Data_NRMSE))
		print 'The minimum NRMSE value found was: '+str(np.min(self.Data_NRMSE))
		print 'The average NRMSE value was: '+str(np.average(self.Data_NRMSE))
		print 'The NRMSE standard deviation was: '+str(np.std(self.Data_NRMSE, dtype=np.float64))
		plt.figure()
		n, bins, patches = plt.hist(self.Data_NRMSE, bins = 'auto', facecolor='blue', alpha=0.75)
		plt.xlabel('NRMSE')
		plt.ylabel('number of files')
		plt.savefig(self.outfile_bestfit_N.replace('.csv','NRMSE_histogram.pdf'))
		plt.close()

	def return_subset(self,electrode='',tissu='',out_name=''):
		if out_name == '':
			subset = ExVivo_Experiment(self.Path,self.Folders,out_name='zmodel_subset_out')
		else:
			subset = ExVivo_Experiment(self.Path,self.Folders,out_name=out_name)
		# look for corresponding data

		# done
		return subset
