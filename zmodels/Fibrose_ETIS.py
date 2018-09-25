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
# First added:  2018-03-20															#
# Last changed: 2018-06-12														#
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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import sys
import os
import glob
import math
import pandas as pd

def xls2csv(folder,xls_radical,increment =0):
	excel_file1 = folder + 'LargeModulo_' + xls_radical + '.xls'
	excel_file2 = folder + 'MediumModulo_' + xls_radical + '.xls'
	excel_file3 = folder + 'LargePhase_' + xls_radical + '.xls'
	excel_file4 = folder + 'MediumPhase_' + xls_radical + '.xls'

	highf_Mod = pd.read_excel(excel_file1, sheet_name='Sheet 1',header = None)
	lowf_Mod = pd.read_excel(excel_file2, sheet_name='Sheet 1',header = None)
	highf_Ph = pd.read_excel(excel_file3, sheet_name='Sheet 1',header = None)
	lowf_Ph = pd.read_excel(excel_file4, sheet_name='Sheet 1',header = None)

	electrodes = ['0','1','2','3','4','5','6','7']

	for electrode in xrange(len(electrodes)):
		outfile = open(folder+'electrode_'+electrodes[electrode]+'_t_'+str( int(xls_radical)+increment).zfill(3) +'.csv','w')
		for k in xrange(5,len(lowf_Mod[0].values)):
			 line =  str(lowf_Mod[0].values[k]) + ',' + str(lowf_Mod[electrode+1].values[k]) + ',' + str(lowf_Ph[electrode+1].values[k]) + '\n'
			 outfile.write(line)
		for k in xrange(16,len(highf_Mod[0].values)):
			 line =  str(highf_Mod[0].values[k]) + ',' + str(highf_Mod[electrode+1].values[k]) + ',' + str(highf_Ph[electrode+1].values[k]) + '\n'
			 outfile.write(line)
		outfile.close()

class ETIS_Experiment(object):
	def __init__(self, Path, out_name='zmodels_out'):
		super(ETIS_Experiment, self).__init__()
		self.Path = Path
		self.out_name = out_name
		self.Data_list = []

		self.files_to_load = []
		self.sorted_files_to_load = []
		# find all csv files
		for infile in glob.glob(os.path.join(self.Path, '*.csv')):
			self.files_to_load.append(infile)
		self.sorted_files_to_load = sorted(self.files_to_load)
		# alpha-sort them to retrieve chronology 
		for in_file_name in self.sorted_files_to_load:
			self.Data_list.append(expe.MeasSignal(in_file_name,'freq'))
		self.Data_NRMSE = np.zeros(len(self.Data_list))

	def plot_all_measurements(self,pict_format ='.pdf'):
		for data in self.Data_list:
			print "saving %s impedance curves"%data.source
			plt.figure()
			plt.loglog(data.freq, data.zimp)
			plt.grid(True, which="both")
			plt.xlabel('frequency (Hz)')
			plt.ylabel('Impedance magnitude')
			plt.savefig(data.source.replace('.csv','_mod'+pict_format))
			plt.close()

			plt.figure()
			plt.semilogx(data.freq, data.zphase)
			plt.grid(True, which="both")
			plt.xlabel('frequency (Hz)')
			plt.ylabel('Impedance magnitude')
			plt.savefig(data.source.replace('.csv','_phase'+pict_format))
			plt.close()

	def identify_models(self,folder=''):
		# files to save results
		self.outfile_RC_N = folder+self.out_name+'_Serial_RC.csv'
		self.outfile_RCPE_N = folder+self.out_name+'_Serial_RCPE.csv'
		self.outfile_bestfit_N = folder+self.out_name+'_bestfit.csv'
		self.folder = folder
		
		# create or empty file, write headers 
		self.outfile_RC_F = open(self.outfile_RC_N,'a')
		self.outfile_RCPE_F = open(self.outfile_RCPE_N,'a')
		self.outfile_bestfit_F = open(self.outfile_bestfit_N,'a')
		self.outfile_RC_F.write('source\t K\t wc \n')
		self.outfile_RCPE_F.write('source\t gamma\tK\t wc  \n')
		self.outfile_bestfit_F.write('source\t model\n')

		for data in self.Data_list:
			print "... identifying on file %s"%data.source
			try:
				# Simple one zero model for initial conditions
				data.least_square_fit('Serial_RC',bounds=([0.0, 0.0],[np.inf, np.inf]))
				# dynamic initial guess for RCPE model
				gamma_dyn = 1
				K_dyn = data.param[-1][0]
				wb_dyn = data.param[-1][1]
				self.outfile_RC_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\n')
				try:
					data.least_square_fit('Serial_RCPE', p0=(gamma_dyn, K_dyn, wb_dyn))
					self.outfile_RCPE_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\n')
				except RuntimeError:
					print 'identification not successful using Least Square algorithm for simple RC model'
					print '... relaunching with gradient descend'
					try:
						data.minimize_error_fit('Serial_RCPE', p0=(gamma_dyn, K_dyn, wb_dyn))
						self.outfile_RCPE_F.write(data.source+'\t'+str(data.param[-1][0])+'\t'+str(data.param[-1][1])+'\t'+str(data.param[-1][2])+'\n')
					except RuntimeError:
						print '... failed to identify Serial_RCPE model'
						self.outfile_RCPE_F.write(data.source+'\t 0.0 \t 0.0 \t 0.0 \n')
			except RuntimeError:
				print 'identification not successful using Least Square algorithm for simple RC model'
				self.outfile_RC_F.write(data.source+'\t 0.0 \t 0.0 \n')
				pass

		# check for the model with the best fit
			models_nrmse = np.zeros(len(data.zimp_ident))
			for k in xrange(len(data.zimp_ident)):
				models_nrmse[k] = nrm.NRMSE(data.zimp,data.zimp_ident[k])
			minNRMSEarg = np.argmin(models_nrmse)
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
		print 'NRMSE results for the set of data used: '
		print 'The maximum NRMSE value found was: '+str(np.max(self.Data_NRMSE))
		print 'The minimum NRMSE value found was: '+str(np.min(self.Data_NRMSE))
		print 'The average NRMSE value was: '+str(np.average(self.Data_NRMSE))
		print 'The NRMSE standard deviation was: '+str(np.std(self.Data_NRMSE, dtype=np.float64))
		self.outfile_RC_F.close()
		self.outfile_RCPE_F.close()
		self.outfile_bestfit_F.close()

	def time_domain_analysis(self, time, unit='hours'):
		self.time = time
		self.time_unit = unit
		self.freq_vector = self.Data_list[0].freq

		# meshgrids
		self.zdB = np.zeros((len(self.time),len(self.freq_vector)))
		self.z = np.zeros((len(self.time),len(self.freq_vector)))
		self.Delta_z = np.zeros((len(self.time),len(self.freq_vector)))
		self.Delta_zdB = np.zeros((len(self.time),len(self.freq_vector)))
		for i in xrange(len(self.time)):
			for j in xrange(len(self.freq_vector)):
				self.zdB[i,j]=self.Data_list[i].zimp_dB[j]
				self.z[i,j]=self.Data_list[i].zimp[j]
				self.Delta_zdB[i,j]=self.Data_list[i].zimp_dB[j]-self.Data_list[0].zimp_dB[j]
				self.Delta_z[i,j]=self.Data_list[i].zimp[j]-self.Data_list[0].zimp[j]
		self.X, self.Y = np.meshgrid(self.time, self.freq_vector)
		# plot meshfrids
		print self.Path+'Zimp_colormap.pdf' 
		fig2 = plot_time_pcolor(self.X, np.log10(self.Y), self.zdB.T,xl='time (hour)',yl='frequency (log10 of Hz)',title='impedance magnitude (dB Ohms)',path=self.Path+'00_Zimp_colormap.pdf')
		print self.Path+'Delta_Zimp_colormap.pdf' 
		fig2 = plot_time_pcolor(self.X, np.log10(self.Y), self.Delta_zdB.T,xl='time (hour)',yl='frequency (log10 of Hz)',title='impedance magnitude difference (dB Ohms)',path=self.Path+'00_Delta_Zimp_colormap.pdf')
		# Evolution of the R-CPE model, compute and save to a csv for later use
		self.R_vector = np.zeros(len(self.Data_list))
		self.OmegaC_vector = np.zeros(len(self.Data_list))
		self.CPEorder_vector = np.zeros(len(self.Data_list))
		self.Q_vector = np.zeros(len(self.Data_list))
		line = 'hours, R, Q, gamma, FC'
		saveparam = open(self.Path+'00_parameters_over_time.csv','w')
		saveparam.write(line)
		for k in xrange(len(self.Data_list)):
			self.CPEorder_vector[k] = self.Data_list[k].param[-1][0]
			self.R_vector[k] = self.Data_list[k].param[-1][1]
			self.OmegaC_vector[k] = self.Data_list[k].param[-1][2]
			self.Q_vector[k] = 1./(self.R_vector[k]*(self.OmegaC_vector[k]**self.CPEorder_vector[k]))
			line = str(self.time[k]) + ',' + str(self.R_vector[k]) + ',' + str(self.Q_vector[k]) + ',' + str(self.CPEorder_vector[k]) + ',' + str(self.OmegaC_vector[k]) + '\n' 
			saveparam.write(line)
		saveparam.close()
		# plot evolution of parameter over time
		fig3 = plt.figure()
		plt.subplot(3, 1, 1)
		plt.plot(self.time, self.R_vector)
		plt.ylabel('R (Ohm)')
		plt.subplot(3, 1, 2)
		plt.plot(self.time, self.OmegaC_vector)
		plt.ylabel('corner frequency (Hz)')
		plt.subplot(3, 1, 3)
		plt.plot(self.time, self.CPEorder_vector)
		plt.ylabel('CPE order')
		plt.savefig(self.Path+'00_parameters_vs_time.pdf')
		# plot evolution of parameter over time
		fig4 = plt.figure()
		plt.subplot(3, 1, 1)
		plt.plot(self.time, self.R_vector)
		plt.ylabel('R (Ohm)')
		plt.subplot(3, 1, 2)
		plt.plot(self.time, self.Q_vector)
		plt.ylabel('Constant Phase Element (F)')
		plt.subplot(3, 1, 3)
		plt.plot(self.time, self.CPEorder_vector)
		plt.ylabel('CPE order')
		plt.savefig(self.Path+'00_physical_parameters_vs_time.pdf')
		# SHOW ALL FIGURES (keep at the end of method)
		#plt.show()

def plot_time_surf(X,Y,Z,xl='',yl='',zl='',path=''):
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,linewidth=1, antialiased=True)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	ax.set_xlabel(xl)
	ax.set_ylabel(yl)
	ax.set_zlabel(zl)
	if not(path==''):
		plt.savefig(path)
	else:
		pass
	return fig, ax

def plot_time_pcolor(X,Y,Z,xl='',yl='',title='',path=''):
	fig = plt.figure()
	surf = plt.pcolormesh(X, Y, Z, cmap=cm.jet,antialiaseds =True,shading='flat')
	cbar = fig.colorbar(surf, shrink=0.5, aspect=5)#, cbarlabel="dB Ohm"
	cbar.set_label("dB Ohm")
	plt.ylabel(yl)
	plt.xlabel(xl)
	plt.title(title)
	if not(path==''):
		plt.savefig(path)
	else:
		pass
	return fig	


