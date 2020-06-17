#!/usr/bin/env python3
# -*- coding: utf -8 -*-
''' -  --------------------------------------------------------------------------------------------
	Program developed, designed and written by: Piotr Zarzycki zarzycki.piotrek@gmail.com
	December 2016 - January 2017
	
	
	
	------------------------------------------------------------------------------------------------
	
	TO DO 
	1) changing number of peaks after plotting - remove/update plots (now the old peaks remains), remove and add simultaneously all figures 
	2) change the increment in the table 
	3) fitting 
	4) physical constraints
	5) remove the plot style options 
	6) keep the number colored in component list 
	7) add physical constrain settings:
		peaks weights  (radioButton) positive,  negative,  all
		peaks position (within data range), everywherer (means  DX +/- 0.5*DX wher DX is range of data in file)
		self-consitent solution (checkButton  = solver in loop will generate various soluton untill it stay unchnaged) 
	
	
'''	

import matplotlib
matplotlib.use('TkAgg')
import numpy  as np
from matplotlib.figure import Figure
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg  
# this works in 2017, now NavigationToolbar2TkAgg class is deprecated 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk # this works


from matplotlib.backend_bases import key_press_handler

import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename
from tkinter import Menu
from scipy.optimize import leastsq
import time, random, sys, getopt
from PIL import Image, ImageTk
from tkinter import messagebox
from tkinter import ttk, BOTH, LEFT, RIGHT, X, Y, TOP

GLOBAL_max_functions    = 20  
GLOBAL_integrator_steps = 1000


def FunctionSet(P, x, nfun, npar, action,function_types): 
	value = 0.0 		# value to return 
	number_peak = 0		# counter of peaks 
	# P 				array of parameters  (expected size nfun * npar)
	# x 				argument value y = f(x)
	# nfun 				number of peaks (functions to fit)
	# npar 				number of parameters in fitting functions (default 3, but code can be modified for arbitrary number)
	# action 			list of options (0,1),  0  - optimize (include function) and 1 - fix (not include function)
	# function_types 	list of string describing function (Gauss, Hubbert, Lorentz)
	#print('inside-ftypes: ',function_types)
	for i in range(0, nfun * npar, npar):
		#print(' in loo fun>',function_types[number_peak],'<')
		if  action[number_peak]==0:
			if  function_types[number_peak][0] == 'G': # check only the first letter 
				#print(' here I am Gauss')
				value += (P[i] * np.exp(-((x-P[i+1])**2)/(2.0*P[i+2])))
			elif function_types[number_peak][0]=='H':
				value += (4.0*P[i] * np.exp(-(x-P[i+1])/P[i+2]) / ((1.0 + np.exp(-(x-P[i+1])/P[i+2]))**2))
			elif function_types[number_peak][0]=='L': # //Lorentz 
				value += (P[i]*( 1.0/ ( 1+ ((x-P[i+1])/P[i+2])**2 )))
			else: # does not recognize function type 
				print('does not recognize the function type ')
		number_peak+=1
	return value
















class Application(tk.Frame):
	def DebugReadPrintUpdateCurrent(self):
		self.ReadTable()
		print('Current table:')
		for i in range(self.number_of_peaks):
			print('f({:2d} {} flag={} par='.format(i+1,self.Current_Function_Types[i], self.Current_Function_Flags[i]),end='')
			for j in range(self.number_of_parameters):
				print(' {:8.5f} {:8.5f} {:8.5f} '.format(self.Current_Parameters[i*3],self.Current_Parameters[i*3+1]\
				,self.Current_Parameters[i*3+2]))
	
	def DebugPrintCurrent(self):
		print('Current.parm:')
		for i in range(self.number_of_peaks):
			print('f({:2d} {} flag={} par='.format(i+1,self.Current_Function_Types[i], self.Current_Function_Flags[i]),end='')
			print(' {:8.5f} {:8.5f} {:8.5f} '.format(self.Current_Parameters[i*3],self.Current_Parameters[i*3+1]\
			,self.Current_Parameters[i*3+2]))
				
			
	def ReadTable(self):
		if len(self.WidgetFixOptChecks) > 0: # check if the table is not empty 
			self.Current_Function_Flags = []
			self.Current_Function_Types = []
			self.Current_Parameters =[] 
			for i in range(self.number_of_peaks):
				self.Current_Function_Types.append(  self.WidgetFunType_Var[i].get() 				)
				self.Current_Function_Flags.append(  self.WidgetFixOptChecks_Var[i].get()			)
				for j in range(self.number_of_parameters):
					self.Current_Parameters.append( self.WidgetParameters_Var[i*self.number_of_parameters +j].get())
				
	def ConvertColorList(self): # convert matplotlib color symbols to tkinter color names 
		for i in range(len(self.PlotsColorList)):
			color_letter = self.PlotsColorList[i]
			if   color_letter == 'b': self.PlotsColorList[i]='blue'
			elif color_letter == 'g': self.PlotsColorList[i]='green'
			elif color_letter == 'c': self.PlotsColorList[i]='cyan'
			elif color_letter == 'r': self.PlotsColorList[i]='red'
			elif color_letter == 'k': self.PlotsColorList[i]='black'
			elif color_letter == 'y': self.PlotsColorList[i]='yellow'
			elif color_letter == 'm': self.PlotsColorList[i]='magenta'
			else: self.PlotsColorList[i]='black'
				 
	
	
	def plotting_table(self): # plot table and also plot data 
		# first remove all 
		for i in range( GLOBAL_max_functions):
			self.PlotsFunctionList[i].set_xdata( [] )
			self.PlotsFunctionList[i].set_ydata( [] )
			self.PlotsFunctionList[i].set_label('')
			self.PlotsFunctionList[i].set_visible(False)
	
		if len(self.xdata) > 0: self.plotting_data() # plot data, flag plot/hide checked inside function 
				
		# read parameters from the Table current - is the parameter table 
		self.ReadTable()
		tmpXX = [xx for xx in self.xdata] # set the x-value span 
		for i in range( len(self.Current_Function_Flags)):
			tmplocal =[ self.Current_Parameters[i*self.number_of_parameters],self.Current_Parameters[i*self.number_of_parameters+1],\
			self.Current_Parameters[i*self.number_of_parameters+2] ]
			tmpOpt = [0] * self.number_of_peaks
			tmpYY = [FunctionSet(tmplocal, xx, 1, 3, tmpOpt ,[self.Current_Function_Types[i],'Unknown']) for xx in tmpXX]
			tmpTotal = [FunctionSet(self.Current_Parameters, xx, self.number_of_peaks, self.number_of_parameters, tmpOpt,\
			self.Current_Function_Types) for xx in  tmpXX]
			if self.CheckShowTotalFit_Var.get() == True:
				self.lineT.set_xdata( tmpXX   )
				self.lineT.set_ydata( tmpTotal )
			else:
				self.lineT.set_xdata( []   )
				self.lineT.set_ydata( []   )
			if self.WidgetPlotChecks_Var[i].get() == True:		
				self.PlotsFunctionList[i].set_xdata( tmpXX )
				self.PlotsFunctionList[i].set_ydata( tmpYY )	
				self.PlotsFunctionList[i].set_label('f# {:2d}'.format(i+1))
				self.PlotsFunctionList[i].set_visible(True)
				self.PlotsFunctionList[i].set_alpha(0.6)
			else:
				self.PlotsFunctionList[i].set_xdata( [] )
				self.PlotsFunctionList[i].set_ydata( [])	
				self.PlotsFunctionList[i].set_label('')
				self.PlotsFunctionList[i].set_visible(False)
				
		self.WidgetLabel[i].configure(foreground= self.PlotsColorList[i])
		self.ax.legend(loc='upper left', bbox_to_anchor=(.0, 1), ncol=1, fancybox=True, shadow=True,fontsize='small')	
		self.ax.axis([self.plot_range[0],self.plot_range[1],self.plot_range[2],self.plot_range[3]])
		self.canvas.draw()	
	
	def calculate_area(self, f, parameters, nfun, npar, action, function_types):
		nsteps = GLOBAL_integrator_steps		# (warning) fixed accuracy of the trapezoidal integrator 
		h_step = (self.data_range[1] - self.data_range[0]) / ( 1.0*nsteps)
		area = (f(parameters,self.data_range[0],nfun,npar,action,function_types) + \
		f(parameters,self.data_range[1],nfun,npar,action, function_types))/2.0
		action=[0]*self.number_of_peaks
		for i in range(nsteps):
			xx = self.data_range[0] + i * h_step
			area += f(parameters, xx, nfun, npar, action, function_types)
		return area * h_step
		
	def CalculateAreasFit(self, parameters, types):
		self.total_surface_area = 0
		self.individual_surface_area = []
		flags_for_area = [0]*self.number_of_peaks 
		self.total_surface_area = self.calculate_area(FunctionSet, parameters, self.number_of_peaks , self.number_of_parameters,\
		 flags_for_area, types)
		for i in range(self.number_of_peaks):
			Pset = [ parameters[i*3],parameters[i*3+1],parameters[i*3+2]]
			self.individual_surface_area.append( self.calculate_area(FunctionSet, Pset, 1, self.number_of_parameters,flags_for_area, [types[i]]*2))
		
	
	def CalculateFWHM(self, parameters, nfun, npar, types): # numerically approximate full width at half maximum (FWHM)
		self.component_FWHM = []
		for i in range(self.number_of_peaks):
			Pset = [ parameters[i*3],parameters[i*3+1],parameters[i*3+2]]
			height = FunctionSet(Pset, Pset[1], 1, self.number_of_parameters, [0]*2, [types[i]]*2)
			nsteps = GLOBAL_integrator_steps # arbitrary accuracy, incrase if needed 
			x  = self.data_range[0]
			dx = (self.data_range[1]-self.data_range[0])/(1.0*nsteps)
			left_x,right_x 		= 0,0
			left_dy, right_dy 	= 1,1
			for j in range(nsteps):
				dy = abs(FunctionSet(Pset, x, 1, self.number_of_parameters, [0]*2, [types[i]]*2) - 0.5*height)
				if(    dy < left_dy  and x < Pset[1]):
					left_x = x; left_dy  = dy
				elif(  dy < right_dy and x > Pset[1]):
					right_x = x; right_dy = dy
				x+=dx
			self.component_FWHM.append( (right_x - left_x) )
		
		
		
	def UpdateTableWithAreasFWHM(self):
		for i in range(self.number_of_peaks):
			self.WidgetComponentFWHM[i].configure(text='{:6.3f}'.format(self.component_FWHM[i]),foreground=self.PlotsColorList[i])
			self.WidgetComponentArea   [i].configure(text='{:10.4f}'.format( self.individual_surface_area[i]),justify=tk.RIGHT)
			self.WidgetComponentWeight [i].configure(justify=tk.RIGHT,\
			text='{:7.3f}'.format( 100.0*self.individual_surface_area[i]/(1.0*self.total_surface_area)))
			self.WidgetComponentArea   [i].configure(foreground=self.PlotsColorList[i])
			self.WidgetComponentWeight [i].configure(foreground=self.PlotsColorList[i])
			self.WidgetLabel[i].configure(foreground=self.PlotsColorList[i])
	
	
	def CopyParameterValues(self, P_from, P_to):
		try:
			for i in range(self.number_of_peaks):
				for j in range(self.number_of_parameters):
					P_to[i*3+j] = P_from[i*3+j]
		except:
			print(' cannot copy parameters ...')
	
	
	
	
	def FitModel(self):
		
		# disable action buttons for the time period of optimization 
		self.FitButton.configure(state='disable')
		self.PlotButton.configure(state='disable')
		self.ResetButton.configure(state='disable')
		
		
		self.final_label.configure(text='')
		self.labelsBottom.configure(text=' Component functions:') 
		
		if self.CheckResetEachTimeVar.get() == True: # check if the existing parameters and settings are used or ignored 
			self.ResetStartingPoint()
			
		random.seed() # update ranodm seed 
		# update ydata with ydata_original
		for i in range(len(self.xdata)):  self.ydata[i] = self.ydata_org[i]
		self.ReadTable()
		peaks_to_fit = self.number_of_peaks
		npeak =0
		for i in range(0,len(self.Current_Parameters),self.number_of_parameters):
			if( self.Current_Function_Flags[npeak]==1): 	# fixed peak 
				tmpParm = [self.Current_Parameters[i],self.Current_Parameters[i+1],self.Current_Parameters[i+2]]
				#print(tmpParm)
				peaks_to_fit-=1
				for i in range(len(self.xdata)):
					self.ydata[i] -= FunctionSet(tmpParm, self.xdata[i], 1, self.number_of_parameters, [0]*2, [self.Current_Function_Types[npeak]]*2)
			npeak+=1
		#print(' peaks = {}, to fit {}'.format( self.number_of_peaks, peaks_to_fit ))
		self.lineT.set_xdata( self.xdata  )
		self.lineT.set_ydata( self.ydata  )
		self.canvas.draw()
		
		self.ReadTable() # write parameters from the table to self.Current_* (Flags, Types, Parameters)
		ErrorFunction=lambda P,x,y: FunctionSet(P,x,self.number_of_peaks,self.number_of_parameters,self.Current_Function_Flags,\
		self.Current_Function_Types)-y
		np.seterr(all='ignore')	# ignore warnings and error messages from scipy 
		
		
		
		# initialize the container to keep the best parameters globally
		theBestParameters = []
		for i in range(self.number_of_peaks):
			for j in range(self.number_of_parameters):
				theBestParameters.append(0.0)
		theBestAccuracy = 1000.0	# initialize the best accuracy with large number		
		
		solver_iteration    =0
		iteration_physical  =False
		solver_converged    =False
		solver_completed    =False
		
		self.ReadTable() # write current table to self.Current_*
		self.CopyParameterValues(self.Current_Parameters, theBestParameters)
		#print(' starting parameters in table')
		#print(theBestParameters,'\n\n')
		
		
		self.progress["value"] = 0
		self.progress["maximum"] = self.max_iterations
		
		self.OptimizationHistory = [] #np.zeros((self.max_iterations, self.number_of_parameters*self.number_of_peaks))
		#print(self.OptimizationHistory)
		how_many_iterations_saved_in_history = 0 
		
		index_of_best_iteration = 0
		
		while ( solver_iteration < (self.max_iterations)  and  solver_completed ==False):
			self.progress_variable.set(solver_iteration+1) 
			#self.label_pogress.configure(text='Progress: {:>4d} Iteration:'.format(solver_iteration),foreground ='steel blue')
			self.label_pogress.configure(text='Iteration # {:6d}:'.format(solver_iteration))#,foreground ='steel blue')
			
			self.update_idletasks()

			solver_iteration+=1
			final_parameters, pcov, infodict, errmsg, success = \
			leastsq(ErrorFunction, self.Current_Parameters, args=(self.xdata, self.ydata), \
			full_output=1, epsfcn=0.0001)
	
			sum_residues_partial =0
			sum_residues_global  =0
			for i in range(len(self.xdata)):
				sum_residues_partial += (abs(self.ydata[i] - FunctionSet(final_parameters, self.xdata[i], self.number_of_peaks, \
				self.number_of_parameters, self.Current_Function_Flags, self.Current_Function_Types))) 
				sum_residues_global  += (abs(self.ydata_org[i] - FunctionSet(final_parameters, self.xdata[i], self.number_of_peaks, \
				self.number_of_parameters, [0]*self.number_of_peaks, self.Current_Function_Types)))
			
			current_fit_accuracy_partial= sum_residues_partial/(1.0*len(self.xdata))
			current_fit_accuracy_global  =  sum_residues_global /(1.0*len(self.xdata))
			self.CopyParameterValues(final_parameters, self.Current_Parameters) 
			
			if( current_fit_accuracy_global <= self.fitting_accuracy  and current_fit_accuracy_partial <= self.fitting_accuracy):
				solver_converged = True 
			else: solver_converged = False 
			iteration_physical = True  # add conditions later
			for i in range(self.number_of_peaks):
				if iteration_physical == True:
					if self.PeakValue.get() ==1 :
						if( self.Current_Parameters[i*3] <0): iteration_physical= False 
					elif self.PeakValue.get() ==2 :
						if( self.Current_Parameters[i*3] >0): iteration_physical= False 
				if  self.PeakHeight.get() ==1 :
					if( self.Current_Parameters[i*3] > self.data_range[3]): iteration_physical= False 
				if  self.PeakPosition.get() ==1 :
					if( self.Current_Parameters[i*3+1] > self.data_range[1] or self.Current_Parameters[i*3+1] < self.data_range[0] ):
						iteration_physical= False
					
			if iteration_physical == True:
				if( current_fit_accuracy_global < theBestAccuracy): # update the best parameters and corresponding accuracy 
					theBestAccuracy = current_fit_accuracy_global
					self.CopyParameterValues(self.Current_Parameters, theBestParameters)
					index_of_best_iteration = solver_iteration
					self.UpdateParameterTable(self.Current_Function_Types, self.Current_Function_Flags, self.Current_Parameters)
					if self.CheckPlotIteration.get() == True : self.plotting_table()
					self.CalculateAreasFit( self.Current_Parameters, self.Current_Function_Types)
					self.CalculateFWHM(self.Current_Parameters, self.number_of_peaks, self.number_of_parameters, self.Current_Function_Types)
					self.UpdateTableWithAreasFWHM()	
					
					'''self.History_Parameters 	 =[] #np.zeros((2,2)) # // meaningless initialization 
					self.History_Accuracy        =[]
					self.History_Total_Area      =[]
					self.History_Component_Area  =[] #np.zeros((2,2)) # // meaningless initialization 
					self.History_theBest         =0	 # index of the best fitting parameters 
					self.History_number_solutions=1  # number of solutions in history 	
					for ip in range(self.number_of_peaks):
						for jp in range(self.number_of_parameters):
							self.OptimizationHistory[ how_many_iterations_saved_in_history ][ip*self.number_of_parameters + jp]=\
							self.Current_Parameters[ip*self.number_of_parameters + jp]
					how_many_iterations_saved_in_history+=1 '''
					
					
			if( solver_iteration % self.shake_iterations == 0 and solver_completed ==False): 
				self.shake_parameters(self.Current_Function_Flags) # shaking the initial parameters 
				self.CopyParameterValues(self.InitialGuess, self.Current_Parameters)	# and copy the shaked parameters to the solverinput Current
				
			if iteration_physical ==False:
				#self.DebugPrintCurrent()
				self.shake_parameters(self.Current_Function_Flags) # shaking the initial parameters 
				self.CopyParameterValues(self.InitialGuess, self.Current_Parameters)	# and copy the shaked parameters to the solverinput Current
				
				
			if( solver_converged == True and iteration_physical == True): 
				solver_completed = True
			else:
				if( solver_iteration % self.restart_iterations ==0 ):
					self.restart_parameters(self.Current_Function_Flags) # again restart initial 
					#print(' restart parameters (more bruttal shake) after ', solver_iteration)
					self.CopyParameterValues(self.InitialGuess, self.Current_Parameters)
					#self.UpdateParameterTable(self.Current_Function_Types, self.Current_Function_Flags, self.InitialGuess) 	
			
		#print( 'how many iterations saved in history = ', how_many_iterations_saved_in_history )
		#for hi in range (how_many_iterations_saved_in_history):
		#	print( ' history iteration ', hi)
		#	print( ' par', self.OptimizationHistory[hi])
			
			#for hf in range (self.number_of_peaks): 
			#	for hp in range (self.number_of_parameters):
			#		print( )
					
		# get final solution and put into table 
		#print(theBestParameters)
		for i in range( self.number_of_peaks):
			for j in range( self.number_of_parameters):
				self.Current_Parameters[i*3+j] = theBestParameters[i*3+j]
		self.progress_variable.set(self.max_iterations) 
		self.update_idletasks()
			
		self.UpdateParameterTable(self.Current_Function_Types, self.Current_Function_Flags, self.Current_Parameters)	
		self.plotting_table()
		
		self.CalculateAreasFit( self.Current_Parameters, self.Current_Function_Types)
		self.CalculateFWHM(self.Current_Parameters, self.number_of_peaks, self.number_of_parameters, self.Current_Function_Types)
		self.UpdateTableWithAreasFWHM()	

		
		if solver_iteration == self.max_iterations:
			self.labelsBottom.configure(text=' Component functions: max. of iteration reached (accuracy= {:6.3e}, requested= {:6.3e})'\
			.format(theBestAccuracy,self.fitting_accuracy ))
		else:		
			self.labelsBottom.configure(text=' Component functions: the best fit among {:3d} iterations, (accuracy {:6.3e} requested {:6.3e})'\
			.format(solver_iteration, theBestAccuracy, self.fitting_accuracy))
		
		self.final_label.configure(text='Best Fit: (iteration {:3d}) out of {:3d} (max={:3d}) with accuracy {:6.3e}'.format\
		(index_of_best_iteration, solver_iteration, self.max_iterations,  theBestAccuracy))
		
		self.FitButton.configure(state='normal')
		self.PlotButton.configure(state='normal')
		self.ResetButton.configure(state='normal')
		
		
	def callback(self,event): # get position of mouse over the point on the graph 
		try: 
			self.POINT_X_LABEL.configure( text='(X={:8.3f}'.format(event.xdata))  
			self.POINT_Y_LABEL.configure( text= 'Y={:8.3f})'.format(event.xdata)) 
		except:
			self.POINT_X_LABEL.configure( text='')  
			self.POINT_Y_LABEL.configure( text='') 
				
		
	def ActivateUserRange(self):
		if self.AutoScalePlot_Var.get() == False:  # no autoscale (unchecked)
			for child in self.FigLabelsRightB.winfo_children():
				child.configure(state='normal')
			self.AutoScalePlot_Obj.deselect()
			
		elif self.AutoScalePlot_Var.get() == True: # autoscale, checked 
			for child in self.FigLabelsRightB.winfo_children():
				child.configure(state='disable')
			self.AutoScalePlot_Obj.configure(state='normal')
			self.AutoScalePlot_Obj.select()
			self.AutoscalePlotRange()
			self.plotting_data()
		
	def ShowHide(self,i):
		ShowPlotFlag = self.WidgetPlotChecks_Var[i].get()
		if ShowPlotFlag == True:
			self.PlotsFunctionList[i].set_label('f# {:2d}'.format(i+1))
			self.PlotsFunctionList[i].set_visible(True)
		elif ShowPlotFlag == False:
			self.PlotsFunctionList[i].set_label('')
			self.PlotsFunctionList[i].set_visible(False)
		self.canvas.draw()
		
		
	def ChangeNumberFunctions(self,event):
		if self.number_of_peaks <= self.selectedNumberFun.get() :
			self.number_of_peaks = self.selectedNumberFun.get()
			self.CreateParameterTable()
		elif self.number_of_peaks > self.selectedNumberFun.get() :
			self.number_of_peaks = self.selectedNumberFun.get()
			self.CreateParameterTable()
			
	def UpdateAccuracy(self,v1,v2,v3): # immediate trace of parameter value 
		try: 
			self.fitting_accuracy = self.RequestedAccuracy.get()
		except:
			self.fitting_accuracy= 0.01
			self.RequestedAccuracy.set( self.fitting_accuracy)
			
				
	def UpdateIterations(self, v1,v2,v3): #trace 
		try:
			self.max_iterations =  self.RequestedIterations.get()
			if self.max_iterations <= 0: raise( ValueError )
		except: 
			self.max_iterations = 500
			self.RequestedIterations.set(self.max_iterations)
			
	def UpdateIterationsRestart(self, v1,v2,v3):
		try:
			self.restart_iterations = self.RequestedRestart.get()
			if self.restart_iterations <=0: raise( ValueError )
		except:
			self.restart_iterations = 50
			self.RequestedRestart.set ( self.restart_iterations)
	
	def UpdateIterationsShake(self, v1,v2,v3):
		try:
			self.shake_iterations = self.RequestedShake.get()
			if self.shake_iterations <=0: raise( ValueError )
		except:
			self.shake_iterations = 5
			self.RequestedShake.set ( self.shake_iterations)
			
	def ResetStartingPoint(self): 
		self.GenerateInitialGuess() 
		self.UpdateParameterTable(self.InitialTypes, self.InitialFlags, self.InitialGuess)
		
		
	def	Create_Toolbar_GUI(self, master): # GENERATE THE TOOLBAR 
		self.toolbar  = tk.Frame(master)  # bg=blue	
		self.toolbar.pack( side=TOP, fill=X)
		toolbar_padx=1
		toolbar_pady=1
		
		try:
			self.icon_fit 		= Image.open('fit_icon.png')
			self.picon_fit 		= ImageTk.PhotoImage(self.icon_fit)
			self.icon_reset 	= Image.open('reset_icon.png')
			self.picon_reset	= ImageTk.PhotoImage(self.icon_reset)
			self.icon_plot 		= Image.open('plot_icon.png')
			self.picon_plot 	= ImageTk.PhotoImage(self.icon_plot)
			self.icon_open_data = Image.open('open_data_icon.png')
			self.picon_open_data= ImageTk.PhotoImage(self.icon_open_data)
			self.icon_open_parm = Image.open('open_parm_icon.png')
			self.picon_open_parm= ImageTk.PhotoImage(self.icon_open_parm)
			self.icon_save_data = Image.open('save_data_icon.png')
			self.picon_save_data= ImageTk.PhotoImage(self.icon_save_data)
			self.icon_save_parm = Image.open('save_parm_icon.png')
			self.picon_save_parm= ImageTk.PhotoImage(self.icon_save_parm)
			#self.icon_save_hist = Image.open('save_history_icon.png')
			#self.picon_save_hist= ImageTk.PhotoImage(self.icon_save_hist)
			self.icon_save_fig 	= Image.open('save_figures_icon.png')
			self.picon_save_fig = ImageTk.PhotoImage(self.icon_save_fig)
			# ------ CREATE BUTTONS WITH IMAGES ---------------------------------------------------------------------------
			self.OpenDataButton = tk.Button(self.toolbar, image=self.picon_open_data, command=self.OpenDataFile)
			self.OpenParametersButton = tk.Button(self.toolbar, image=self.picon_open_parm, command=self.OpenParameterFile)
			self.SaveDataButton = tk.Button(self.toolbar, image=self.picon_save_data,command=self.SaveDataFile)
			self.SaveParametersButton = tk.Button(self.toolbar, image=self.picon_save_parm,command=self.SaveParameterFile)
			#self.SaveHistoryButton = tk.Button(self.toolbar,image=self.picon_save_hist)
			self.SaveFigureButton = tk.Button(self.toolbar, image=self.picon_save_fig,command=self.SaveFigureAll)
			self.FitButton = tk.Button(self.toolbar,image=self.picon_fit,command=self.FitModel) 
			self.PlotButton = tk.Button(self.toolbar ,command=self.plotting_table,image=self.picon_plot) 
			self.ResetButton = tk.Button(self.toolbar,command=self.ResetStartingPoint,image=self.picon_reset) 
		except: 
			# ------ CREATE BUTTONS WITH NAMES (IMAGES NOT FOUND) ----------------------------------------------------------
			self.OpenDataButton = tk.Button(self.toolbar, text='Open Data', command=self.OpenDataFile)
			self.OpenParametersButton = tk.Button(self.toolbar, text='Open Parameters', command=self.OpenParameterFile)
			self.SaveDataButton = tk.Button(self.toolbar, text='Save Data',command=self.SaveDataFile)
			self.SaveParametersButton = tk.Button(self.toolbar, text='Save Parameters',command=self.SaveParameterFile)
			#self.SaveHistoryButton = tk.Button(self.toolbar,text='Save History')
			self.SaveFigureButton = tk.Button(self.toolbar, text='Save Figure',command=self.SaveFigureAll)
			self.FitButton = tk.Button(self.toolbar, text='Fit' ,command=self.FitModel) 
			self.PlotButton = tk.Button(self.toolbar ,command=self.plotting_table, text='Plot') 
			self.ResetButton = tk.Button(self.toolbar,command=self.ResetStartingPoint, text='Reset') 
		# ------ ARRANGE BUTTONS -------------------------------------------------------------------------
		self.POINT_X_LABEL = ttk.Label( self.toolbar, text='',font='Helvetica 12')
		self.POINT_Y_LABEL = ttk.Label( self.toolbar, text='',font='Helvetica 12')
		
		
		self.OpenDataButton.pack(side=tk.LEFT,padx=toolbar_padx, pady=toolbar_pady)
		self.OpenParametersButton.pack(side=tk.LEFT,padx=toolbar_padx, pady=toolbar_pady)
		self.SaveDataButton.pack(side=tk.LEFT,padx=toolbar_padx, pady=toolbar_pady)
		self.SaveParametersButton.pack(side=tk.LEFT,padx=toolbar_padx, pady=toolbar_pady)
		#self.SaveHistoryButton.pack(side=tk.LEFT,padx=toolbar_padx, pady=toolbar_pady)
		self.SaveFigureButton.pack(side=tk.LEFT,padx=toolbar_padx, pady=toolbar_pady)
		
		self.FitButton.pack(side=RIGHT,padx=toolbar_padx, pady=toolbar_pady)
		self.PlotButton.pack(side=RIGHT,padx=toolbar_padx, pady=toolbar_pady)
		self.ResetButton.pack(side=RIGHT,padx=toolbar_padx, pady=toolbar_pady)
		
		self.POINT_Y_LABEL.pack(side=RIGHT, padx=12, pady=toolbar_pady)
		self.POINT_X_LABEL.pack(side=RIGHT, padx=12, pady=toolbar_pady)
		
		
		
	def Create_Sides_GUI(self, master):
		# labels and widgets 
		self.WindowLeft = tk.Frame(master)
		self.WindowLeft.pack(side=LEFT,expand=True,fill=BOTH,padx=2,pady=5)
		self.WindowRight = tk.Frame(master)
		self.WindowRight.pack(side=LEFT,expand=True,fill=BOTH,padx=2,pady=5)

		self.labelsTop = ttk.LabelFrame(self.WindowLeft, text='Settings ')
		self.labelsTop.pack(side=TOP,expand=True,fill=BOTH)
		
		self.labelsMiddle = ttk.LabelFrame(self.WindowLeft, text=' Fitting Progress ')
		self.labelsMiddle.pack(side=TOP,expand=True,fill=BOTH,padx=2, pady=5)
		
		self.labelsBottom = ttk.LabelFrame(self.WindowLeft, text=' Components : ')
		self.labelsBottom.pack(side=TOP,expand=True,fill=BOTH,padx=2,pady=2)
		
		ttk.Label(self.labelsTop, text='Function Type:').grid(column=0,row=0,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='Number of Functions:').grid(column=0,row=1,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='Fitting Accuracy:').grid(column=0,row=2,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='Max Iterations:').grid(column=0,row=3,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='Solver: Shake after ',foreground='sea green').grid(column=0,row=4,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='Solver: Restart after ',foreground='sea green').grid(column=0,row=5,sticky=tk.W,padx=5,pady=2)
		
		ttk.Label(self.labelsTop, text='Constraints:',foreground ='steel blue',font='bold')\
		.grid(column=3,row=0,columnspan=4,sticky=tk.N,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='function values:',foreground ='steel blue').grid(column=3,row=1,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='peak position:',foreground ='steel blue').grid(column=3,row=2,sticky=tk.W,padx=5,pady=2)
		ttk.Label(self.labelsTop, text='peak height',foreground ='steel blue').grid(column=3,row=3,sticky=tk.W,padx=5,pady=2)
		self.label_pogress=ttk.Label(self.labelsMiddle, text='Iteration # {:6d}:'.format(0))#,foreground ='steel blue')
		self.label_pogress.grid(column=0,row=0,sticky=tk.W,padx=5,pady=2)
		self.final_label = ttk.Label(self.labelsMiddle, text='')
		self.final_label.grid(column=0,row=1,columnspan=3,sticky=tk.W,padx=5,pady=2)
		
		self.PeakValue = tk.IntVar()
		self.RadPeak_Obj_Positive = tk.Radiobutton (self.labelsTop,text='positive',variable=self.PeakValue,value=1)
		self.RadPeak_Obj_Negative = tk.Radiobutton (self.labelsTop,text='negative',variable=self.PeakValue,value=2)
		self.RadPeak_Obj_Any = tk.Radiobutton (self.labelsTop,text='any',variable=self.PeakValue,value=3)
		self.RadPeak_Obj_Positive.grid(column=4,row=1,sticky=tk.W,padx=0,pady=1)
		self.RadPeak_Obj_Negative.grid(column=5,row=1,sticky=tk.W,padx=0,pady=1)
		self.RadPeak_Obj_Any	 .grid(column=6,row=1,sticky=tk.W,padx=0,pady=1)
		self.PeakValue.set(1)
		
		self.PeakPosition = tk.IntVar()
		self.PeakPosition_within 	= tk.Radiobutton (self.labelsTop,text='inside data',variable=self.PeakPosition,value=1)
		self.PeakPosition_anywhere  = tk.Radiobutton (self.labelsTop,text='any',variable=self.PeakPosition,value=2)
		self.PeakPosition_within	.grid(column=4,row=2,sticky=tk.W,padx=0,pady=1)
		self.PeakPosition_anywhere	.grid(column=5,row=2,sticky=tk.W,padx=0,pady=1)
		self.PeakPosition.set(1)
		
		self.PeakHeight = tk.IntVar()
		self.Peakheight_below  	= tk.Radiobutton (self.labelsTop,text='below max(y)',variable=self.PeakHeight,value=1)
		self.Peakheight_any   = tk.Radiobutton (self.labelsTop,text='any',variable=self.PeakHeight,value=2)
		self.Peakheight_below 	.grid(column=4,row=3,sticky=tk.W,padx=0,pady=1)
		self.Peakheight_any 	.grid(column=5,row=3,sticky=tk.W,padx=0,pady=1)
		self.PeakHeight.set(1)
		
		self.progress_variable = tk.DoubleVar()
		self.progress = ttk.Progressbar(self.labelsMiddle, orient="horizontal", variable= self.progress_variable, length=425)
		self.progress.grid(column=1,row=0,columnspan=4,sticky=tk.W,padx=5,pady=1)
		self.label_pogress.configure(text='Iteration # {:6d}:'.format(0),foreground ='steel blue')
		self.selectFun_ComboBox = ttk.Combobox(self.labelsTop,width=8,textvariable=self.selected_function_type,state='readonly',justify=tk.CENTER)
		self.selectFun_ComboBox['values']=('Gauss','Hubbert','Lorentz')
		self.selectFun_ComboBox.grid(column=1,row=0,sticky=tk.E,padx=5,pady=2)
		self.selectFun_ComboBox.current(0)
		
		self.selectedNumberFun = tk.IntVar()
		self.selectNumberFun_ComboBox = ttk.Combobox(self.labelsTop,width=8,justify=tk.RIGHT,textvariable=self.selectedNumberFun)
		self.selectNumberFun_ComboBox['values']=tuple([ i for i in range(1,GLOBAL_max_functions+1)])
		self.selectNumberFun_ComboBox.grid(column=1,row=1,sticky=tk.E,padx=5,pady=2)
		self.selectNumberFun_ComboBox.current(4)
		
		self.selectNumberFun_ComboBox.bind("<<ComboboxSelected>>", self.ChangeNumberFunctions)
		self.RequestedAccuracy=tk.DoubleVar()
		self.AccuracySpinbox = tk.Spinbox(self.labelsTop, from_ = 0.0001, to = 0.1,  increment=0.0001, width=8, format="%6.5f",justify=tk.RIGHT,\
		textvariable=self.RequestedAccuracy)
		self.AccuracySpinbox.grid(column=1,row=2,sticky=tk.E,padx=5,pady=2)
		self.RequestedAccuracy.set(self.fitting_accuracy)
		self.RequestedAccuracy.trace('w',self.UpdateAccuracy)
		
		self.RequestedIterations=tk.IntVar()
		self.IterationSpin = tk.Spinbox(self.labelsTop, from_ = 1, to = 1000,  increment=1, width=8,justify=tk.RIGHT,\
		textvariable=self.RequestedIterations)
		self.IterationSpin.grid(column=1,row=3,sticky=tk.E,padx=5,pady=2)
		self.RequestedIterations.set(self.max_iterations)
		self.RequestedIterations.trace('w',self.UpdateIterations)
		
		self.RequestedShake=tk.IntVar()
		self.ShakeSpin = tk.Spinbox(self.labelsTop, from_ = 1, to = 1000,  increment=1, width=8,justify=tk.RIGHT,textvariable=self.RequestedShake)
		self.ShakeSpin.configure( foreground='sea green')
		self.ShakeSpin.grid(column=1,row=4,sticky=tk.E,padx=5,pady=2)
		self.RequestedShake.set(self.shake_iterations)
		self.RequestedShake.trace('w',self.UpdateIterationsShake)
		
		self.RequestedRestart=tk.IntVar()
		self.RestartSpin = tk.Spinbox(self.labelsTop,from_=1,to=1000,increment=1,width=8,justify=tk.RIGHT,textvariable=self.RequestedRestart)
		self.RestartSpin.configure( foreground='sea green')
		self.RestartSpin.grid(column=1,row=5,sticky=tk.E,padx=5,pady=2)
		self.RequestedRestart.set(self.restart_iterations)
		self.RequestedRestart.trace('w',self.UpdateIterationsRestart)
		
		self.CheckPlotIteration=tk.BooleanVar()
		self.CheckPlotIteration_Obj = tk.Checkbutton(self.labelsTop,variable=self.CheckPlotIteration,text='plot each iterative solution')	
		self.CheckPlotIteration_Obj.grid(column=3,row=4,columnspan=3,sticky=tk.W,padx=1,pady=2)
		self.CheckPlotIteration_Obj.select()
		
		self.CheckSaveHistory_Obj = tk.Checkbutton(self.labelsTop,variable=self.CheckPlotIteration,text='save optimization history')	
		self.CheckSaveHistory_Obj.grid(column=3,row=5,columnspan=3,sticky=tk.W,padx=1,pady=2)
		
		self.CheckResetEachTimeVar = tk.BooleanVar()
		self.CheckResetEachTime_Obj = tk.Checkbutton(self.labelsTop,variable=self.CheckResetEachTimeVar,\
		text='ignore current solution (reset guess before fitting)')	
		self.CheckResetEachTime_Obj.grid(column=0,row=6,columnspan=6,sticky=tk.W,padx=1,pady=2)
		self.CheckResetEachTime_Obj.deselect()
		
		#self.CheckAutoUpdateAreasTable_Var = tk.BooleanVar()
		#self.CheckAutoUpdateAreasTable_Obj = tk.Checkbutton(self.labelsTop,variable=self.CheckAutoUpdateAreasTable_Var,\
		#text='autoupdate Area,%,FWHM')
		
		self.CheckShowDataPlot_Var = tk.BooleanVar()
		self.CheckShowDataPlot_Obj = tk.Checkbutton(self.labelsTop,variable=self.CheckShowDataPlot_Var,text='plot data',command=self.plotting_data)
		self.CheckShowDataPlot_Obj.grid(column=5,row=5,columnspan=3,sticky=tk.W,padx=1,pady=2)
		self.CheckShowDataPlot_Obj.select()
		
		self.CheckShowTotalFit_Var = tk.BooleanVar()
		self.CheckShowTotalFit_Obj = tk.Checkbutton(self.labelsTop,variable=self.CheckShowTotalFit_Var,text='plot total fit',\
		command=self.plotting_table)
		self.CheckShowTotalFit_Obj.grid(column=5,row=4,columnspan=3,sticky=tk.W,padx=1,pady=2)
		self.CheckShowTotalFit_Obj.select()
		
		# ------------------------------- CREATE FIGURE -----------------------------------------------------------------------	
		self.fig  = Figure(facecolor='white') 
		self.ax   = self.fig.add_subplot(111)
		self.line1, = self.ax.plot([],[],'ko',color='w',label='data')
		self.lineT, = self.ax.plot([],[],'',color='k', lw=3,label='total fit')
		for i in range(GLOBAL_max_functions):
			self.PlotsFunctionList.append('fplot_'	+str(i))
		for i in range( GLOBAL_max_functions ):
			self.PlotsFunctionList[i], = self.ax.plot([],[],lw=3,label='f#{:2d}'.format( i+1),alpha=0.8)
			self.PlotsFunctionList[i].set_visible(False)
			self.PlotsFunctionList[i].set_label('')
			self.PlotsColorList.append( self.PlotsFunctionList[i].get_color() )
		
		self.ConvertColorList()
		self.ax.grid(); 	
		self.ax.legend(loc='upper left', bbox_to_anchor=(.0, 1), ncol=1, fancybox=True, shadow=True,fontsize='small')
		self.canvas = FigureCanvasTkAgg(self.fig,master=self.WindowRight)
		#self.canvas.show()  # this is replaced by .draw() now see line below 
		self.canvas.draw()
		self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
		self.fig.canvas.callbacks.connect('button_press_event', self.callback)
		self.CreateParameterTable()
		
		self.FigLabelsRightB= ttk.LabelFrame(self.WindowRight,text='Plot range:')
		self.FigLabelsRightB.pack(fill=tk.BOTH,expand=True)
		
		ttk.Label(self.FigLabelsRightB, text='Plot ranges: x from ').grid(column=0,row=0,sticky=tk.W,padx=2,pady=5)
		ttk.Label(self.FigLabelsRightB, text=' to ').grid(column=2,row=0,sticky=tk.W,padx=2,pady=5)
		ttk.Label(self.FigLabelsRightB, text=' y from ').grid(column=4,row=0,sticky=tk.W,padx=2,pady=5)
		ttk.Label(self.FigLabelsRightB, text=' to ').grid(column=6,row=0,sticky=tk.W,padx=2,pady=5)
		
		self.AutoScalePlot_Var= tk.BooleanVar()
		self.AutoScalePlot_Obj= tk.Checkbutton(self.FigLabelsRightB,variable=self.AutoScalePlot_Var,text='autoscale',command=self.ActivateUserRange)
		self.AutoScalePlot_Obj.grid(column=8,row=0,sticky=tk.W,padx=1,pady=1)
		
		self.UserPlotRange				 =  []  # values for plot range entered by user 
		self.UserEntryPlotRange_Obj      =  []  # lis of EntryBox for user entry of range (Widgets)
		for j in range(4): self.UserPlotRange.append( tk.DoubleVar()) 
		for j in range(4): self.UserEntryPlotRange_Obj.append(ttk.Entry(self.FigLabelsRightB,textvariable=self.UserPlotRange[j],width=6)) 
		for j in range(4): self.UserEntryPlotRange_Obj[j].grid(column=2*j+1,row=0,sticky=tk.W,padx=1,pady=1)
		for j in range(4): self.UserEntryPlotRange_Obj[j].bind("<Return>",self.UserChangedRange)
		
		for child in self.FigLabelsRightB.winfo_children():
			child.configure(state='disable')
		self.AutoScalePlot_Obj.configure(state='normal')
		self.AutoScalePlot_Obj.select()	
		
			
		
		
		
		
		
		
		
		
	def __init__(self, master=None):
		
		super().__init__(master)
		self.master = master
		
		# GLOBAL VARIABLES FOR MODEL 
		self.FittingParameters      = []
		self.update_parameter_table = False
		self.ShowDataPlot = True
		self.ShowTotalFitPlot = True 


		self.WidgetParameters_Obj	 =[]
		self.WidgetFixOptChecks  	 =[]
		self.WidgetFixOptChecks_Var  =[]
		
		# ---------------------------------------------------------------------------------------
		self.History_Parameters 	 =[] #np.zeros((2,2)) # // meaningless initialization 
		self.History_Accuracy        =[]
		self.History_Total_Area      =[]
		self.History_Component_Area  =[] #np.zeros((2,2)) # // meaningless initialization 
		self.History_theBest         =0	 # index of the best fitting parameters 
		self.History_number_solutions=1  # number of solutions in history 		
		
				
		self.fitting_accuracy   	 = 0.01
		self.max_iterations     	 = 500
		self.shake_iterations   	 =  5
		self.restart_iterations 	 = 50
		
		self.data_range  =[] #min(self.xdata),max(self.xdata),min(self.ydata),max(self.ydata)]
		self.plot_range  =[] # this is used to set axis 
		# remove checks for plots 
		self.WidgetPlotRemChecks  	 	=[]
		self.WidgetPlotRemChecks_Var 	=[]
		self.WidgetComponentArea        =[]
		self.WidgetComponentWeight		=[]
		self.WidgetComponentFWHM		=[]
		
		self.WidgetIncremenet_Var    =[]
		self.WidgetIncrementEntry   =[]		
		self.WidgetParameters_Var	 =[]
		self.WidgetPlotChecks    	 =[]
		self.WidgetPlotChecks_Var	 =[]
		self.WidgetLabel 		 	 =[]
		self.WidgetFunType_Obj      =[]
		self.WidgetFunType_Var      =[]
		self.MinMaxDataValues		 =[] # minX maxX minY maxY of input data (data for fitting) 
		self.pack(fill=BOTH, expand=True)
		self.PlotsFunctionList      =[]
		self.PlotsColorList         =[] # letters representing colors of plot 
		self.InitializeMenu()
		self.xdata = []				# input data  x
		self.ydata = []				# input data  y (this may be changed in partial fit)
		self.ydata_org =[]			# input data y reference (this remains unchanged)
		self.input_data_file =''
		self.flag_plot_data = True
		self.number_of_peaks = 5
		self.number_of_peaks_in_table = self.number_of_peaks
		self.x_at_maximal_y_value = 0 
		self.funtion_types =[] # string in 'Gauss', 'Hubbert', 'Lorentz'
		self.action_type   =[] # 'opt' or 'fix'   (optimize or fixed) 
		self.total_surface_area    = 0		# total surface area under the sum-of-functions (global fit)
		self.component_surface_area = []	# surface area of individual functions 
		self.component_contribution = []	# % contribuion of a given function 
		self.component_FWHM			= []
		self.number_of_parameters = 3 # default value for Guassian peaks 
		self.selected_function_type = tk.StringVar()
		self.result1       = tk.StringVar()
		self.flag_realistic_parameters = True 
		self.ParameterIncrease  = [0.01, 0.1, 0.1]
		self.User_ParmGuess = []
		self.User_Parm_Flag = []
		self.User_FuncType  = []
		self.Current_Parameters     =[] # array of parameter values size = number_of_functions * number_of_parameters
		self.Current_Function_Types =[] # list of function types (Gauss, Hubbert, Lorentz)
		self.Current_Function_Flags =[] # flag of action (fix or opt)
		self.Current_Plot_Flag 	    =[] # flag of action (fix or opt)  # upskirt 
		self.InitialGuess 		= [] # initial guess if none is provided by User 
		self.InitialTypes       = [] # inital types of function 
		self.InitialFlags		= [] # initial flags  
		self.DistributionRange = []  # range of x and y values that contains majority of data 
		self.GuessScale        = []  # scale step for n-initial functions to guess
		self.Create_Toolbar_GUI(master)
		self.Create_Sides_GUI(master)
		
		
		
		
		
			
	def UserChangedRange(self,event):	
		self.UpdatePlotRange()
		self.plotting_data() 	# 
		
	def GenerateInitialGuess(self): # 
		self.InitialFlags =[]
		self.InitialTypes =[]
		self.InitialFlags=[0]*self.number_of_peaks
		self.InitialTypes=[ self.selected_function_type.get() ]*self.number_of_peaks # one type for all 
		self.InitialGuess = [] # reset the initial guess 
		
		if len(self.xdata) >0:
			self.GuessScale   = (self.DistributionRange[2] -self.DistributionRange[0])/(1.0*self.number_of_peaks)
			self.InitialGuess.append( self.data_range[3]  			) 	# max of y 
			self.InitialGuess.append( self.x_at_maximal_y_value    )	# P.append( maxX )
			self.InitialGuess.append( random.uniform( 0.5, 1.0)    )	# random between 0.5 and 1.0 
			for i in range( self.number_of_peaks - 1):	# for the other peaks random position with the spectrum max-value 
				self.InitialGuess.append( self.data_range[3]  			) 	# max of y 
				self.InitialGuess.append( random.uniform( self.DistributionRange[0], self.DistributionRange[2]))
				self.InitialGuess.append( random.uniform( 0.5, 1.0)   )	# random between 0.5 and 1.0 
		else:
			self.InitialGuess =[]
			for i in range( self.number_of_peaks):
					self.InitialGuess.append(random.uniform(0.2,1.0))
					self.InitialGuess.append(random.uniform(1.0,5.0))
					self.InitialGuess.append(random.uniform(0.1,10.))
					
				
	def	shake_parameters(self,action):
		#random.seed()
		for i in range(self.number_of_peaks):
			if action[i]==0:	# check if peak is not fitted 
				self.InitialGuess[i*3  ] = self.data_range[3]  				# max of y 
				self.InitialGuess[i*3+1] = random.uniform( self.DistributionRange[0], self.DistributionRange[2])
				self.InitialGuess[i*3+2] = random.uniform( 0.5, 1.0)   	# random between 0.5 and 1.0
			else:
				self.InitialGuess[i*3  ] = self.Current_Parameters[i*3]  				
				self.InitialGuess[i*3+1] = self.Current_Parameters[i*3+1]
				self.InitialGuess[i*3+2] = self.Current_Parameters[i*3+2]
		
							

	def	restart_parameters(self,action): 
		random.seed()  
	# similar to shake but with peaks anywhere within the data range and less resctrictive with respect to height and peak position 
		for i in range(self.number_of_peaks):
			if action[i]==0:	# check if peak is not fitted 
				self.InitialGuess[i*3  ] = random.uniform(0.1*self.data_range[3], self.data_range[3]) # max of y 
				self.InitialGuess[i*3+1] = random.uniform( self.data_range[0], self.data_range[1])
				self.InitialGuess[i*3+2] = random.uniform( 0.5, 10.0)   		# random between 0.2 and 10.0
			else:
				self.InitialGuess[i*3  ] = self.Current_Parameters[i*3]  				
				self.InitialGuess[i*3+1] = self.Current_Parameters[i*3+1]
				self.InitialGuess[i*3+2] = self.Current_Parameters[i*3+2]
		#self.UpdateParameterTable(self.Current_Function_Types, self.Current_Function_Flags, self.InitialGuess)

	
			

	def AutoscalePlotRange(self): # autoscale the axis on the figure (based on the x,y range of data)
		for j in range(4): 
			self.UserPlotRange[j].set( '{:5.2f}'.format(self.data_range[j]) )
			self.plot_range[j] = self.data_range[j]
		self.plot_range[3] += 0.05*(abs(self.data_range[3]-self.data_range[2])) # increase the y-axis by 5% of spectrum span for visual astetics 
		self.plot_range[2] -= 0.05*(abs(self.data_range[3]-self.data_range[2])) #  
		self.plotting_data()	
		
	
	
	def UpdatePlotRange(self): # change the axis on the figure according to User input 
		for j in range(4): 
			self.plot_range[j]=self.UserPlotRange[j].get()
		
	def SaveDataFile(self):
		try: # save data to file 
			self.file_to_save_data = asksaveasfilename(initialdir="",filetypes =(("Data File","*.dat"),\
			("Text File", "*.txt"),("All Files","*.*")),title = "File name to save data.")
			self.ReadTable()
			print(' filename =',self.file_to_save_data)
			output_file = open(self.file_to_save_data,'w')
			print('#   X(exp)     Y(exp)    F(total)  ',end='',file=output_file)	
			for i in range(self.number_of_peaks):
				print('    f({:2d})  '.format(i+1),end='',file=output_file)
			print('',file=output_file)
			
			tmpOpt = [0]*self.number_of_peaks
			for i in range(len( self.xdata )):
				print('{:10.5f} {:10.5f} {:10.5f} '.format( self.xdata[i], self.ydata_org[i], FunctionSet(self.Current_Parameters, self.xdata[i],\
			 	self.number_of_peaks, self.number_of_parameters, tmpOpt,self.Current_Function_Types)),end='',file=output_file)
				for j in range(self.number_of_peaks):
					tmpPar = [self.Current_Parameters[j*3],self.Current_Parameters[j*3+1],self.Current_Parameters[j*3+2]]
					print(' {:10.5f}'.format(FunctionSet(tmpPar, self.xdata[i], 1,self.number_of_parameters,\
					tmpOpt,self.Current_Function_Types)),end='',file=output_file)
				print('',file=output_file)
			output_file.close()
		except: # 
			pass 
				
		
	def SaveParameterFile(self): # add more details of fitting (accuracy, data)
		try: 
			self.file_to_save_parm = asksaveasfilename(initialdir="", \
			filetypes =(("fitting","*.fit"),("Text File", "*.txt"),("Data File","*.dat"),("All Files","*.*")),title = "File to name.")
			self.ReadTable()	
			ofile_parameters = open(self.file_to_save_parm,'w');
			for i in range(self.number_of_peaks):
				print(" {}  {:5d}  {:8.5f}  {:8.5f}  {:8.5f}  fix  ".\
				format(self.Current_Function_Types[i],i+1,self.Current_Parameters[i*3],self.Current_Parameters[i*3+1],\
				self.Current_Parameters[i*3+2]), file=ofile_parameters)
			ofile_parameters.close()	
		except: 
			pass 
	
	def SaveFigurePDF(self):
		try: 
			self.file_to_save_fig = asksaveasfilename(initialdir="", \
			filetypes =(("PDF File","*.pdf"),("pdfs File", "*.PDF")),title = "File to name.")
			self.fig.savefig(self.file_to_save_fig, dip=600,bbox_inches='tight')
		except: 
			pass
		

		
	def SaveFigurePNG(self):
		try:
			self.file_to_save_fig = asksaveasfilename(initialdir="", \
			filetypes =(("PNG File", "*.png"),("png File","*.PNG")),title = "File to name.")
			self.fig.savefig(self.file_to_save_fig, dip=600,bbox_inches='tight')
		except:
			pass 
		
		
	def SaveFigureSVG(self):
		try:
			self.file_to_save_fig = asksaveasfilename(initialdir="", \
			filetypes =(("SVG File","*.svg"),("Vector Graphics","*.SVG")),title = "File to name.")
			self.fig.savefig(self.file_to_save_fig, dip=600,bbox_inches='tight')
		except:
			pass 
	
		
	def ShowInfo(self):
		toplevel = tk.Toplevel(background="white")
		toplevel.title('About Decomposer (ver. 1.0, January 2017)')
		toplevel.configure(bg='white')
		
		try:
			toplevel.geometry('{}x{}'.format(405, 320))
			image = Image.open("logo_piotrek.png")	# file with the program logo 
			
			photo = ImageTk.PhotoImage(image)
			ttk.Label(toplevel,image=photo,justify=tk.CENTER,background="white").grid(column=0,row=0) # show graphics 
		except: 
			tk.Label(toplevel,text=" Decomposer v1.0  (January, 2017)",justify=tk.CENTER,bg='white').grid(column=0,row=0)
			tk.Label(toplevel,text=" Piotr Zarzycki, Kevin M. Rosso, ...",justify=tk.CENTER,bg='white').grid(column=0,row=1,columnspan=2)
			tk.Label(toplevel,text=" ...", justify=tk.CENTER,bg='white').grid(column=0,row=2,padx=10,pady=10)
			tk.Label(toplevel,text=" ... ",justify=tk.CENTER,bg='white').grid(column=0,row=3,padx=10,pady=10)
			
			
			
			
			pass 
			# labels without image instead of image 
		
		
		toplevel.resizable(0,0)
		toplevel.focus_set()
		toplevel.mainloop()
		
		
		
	def InitializeMenu(self):
		self.master.title("Decomposer 1.0: Spectral Decomposition into bell-shaped functions")
		menuBar = Menu(self.master)
		self.master.config(menu=menuBar)
		
		fileMenu = Menu(menuBar, tearoff=0)
		fileMenu.add_command(label="Open data",command=self.OpenDataFile)
		fileMenu.add_command(label="Open parameters",command=self.OpenParameterFile)
		fileMenu.add_command(label="Save results",command=self.SaveDataFile)
		fileMenu.add_command(label="Save fitting parameters",command=self.SaveParameterFile) 
		fileMenu.add_separator()
		fileMenu.add_command(label="Exit", command=root.destroy)
		fileMenu.add_separator()
		menuBar.add_cascade(label="File", menu=fileMenu)
		
		figureMenu = Menu(menuBar, tearoff=0)
		figureMenu.add_command(label="Save figure to PDF",command=self.SaveFigurePDF)
		figureMenu.add_command(label="Save figure to PNG",command=self.SaveFigurePNG)
		figureMenu.add_command(label="Save figure to SVG",command=self.SaveFigureSVG)
		menuBar.add_cascade(label="Figure", menu=figureMenu)

		actionMenu = Menu(menuBar, tearoff=0)
		actionMenu.add_command(label='Fit data',command=self.FitModel)
		actionMenu.add_command(label='Plot table',command=self.plotting_table)
		actionMenu.add_command(label='Reset initial guess',command=self.ResetStartingPoint)
		menuBar.add_cascade(label='Action',menu=actionMenu)
						
		
		# Add another Menu to the Menu Bar and an item
		helpMenu = Menu(menuBar, tearoff=0)
		helpMenu.add_command(label="Help")
		helpMenu.add_command(label="About",command=self.ShowInfo)
		menuBar.add_cascade(label="Help", menu=helpMenu)	
		
		
	
	def Show_Parameters_GUI(self):
		pass
	
	
	
	def UpdateIncrement(self,v1,v2,v3):	# trace
		try:
			for i in range(self.number_of_parameters):
				self.ParameterIncrease[i] = self.WidgetIncremenet_Var[i].get()
			for i in range(self.number_of_peaks):
				for j in range(self.number_of_parameters):
					self.WidgetParameters_Obj[i*3+j].configure(increment=self.ParameterIncrease[j]) 
		except: 
			for i in range(self.number_of_parameters):
				self.ParameterIncrease[i] = 0.01
				self.WidgetIncremenet_Var[i].set(self.ParameterIncrease[i])
			for i in range(self.number_of_peaks):
				for j in range(self.number_of_parameters):
					self.WidgetParameters_Obj[i*3+j].configure(increment=self.ParameterIncrease[j]) 
		 
			
	
	
	def UpdateParametersUser(self,v1,v2,v3): # trace 
		self.ReadTable()
		self.plotting_table()
		self.CalculateAreasFit( self.Current_Parameters, self.Current_Function_Types)
		self.CalculateFWHM(self.Current_Parameters, self.number_of_peaks, self.number_of_parameters, self.Current_Function_Types)
		self.UpdateTableWithAreasFWHM()
		
	def UpdateParametersUser(self,event): # bind 
		self.ReadTable()
		#print(' total area before = ',self.total_surface_area)
		self.CalculateAreasFit( self.Current_Parameters, self.Current_Function_Types)
		self.CalculateFWHM(self.Current_Parameters, self.number_of_peaks, self.number_of_parameters, self.Current_Function_Types)
		self.UpdateTableWithAreasFWHM()
		#print(' total area after = ',self.total_surface_area)
		self.plotting_table()
			
	
	def plotting_table_and_update(self):
		self.ReadTable()
		#print(' total area before = ',self.total_surface_area)
		self.CalculateAreasFit( self.Current_Parameters, self.Current_Function_Types)
		self.CalculateFWHM(self.Current_Parameters, self.number_of_peaks, self.number_of_parameters, self.Current_Function_Types)
		self.UpdateTableWithAreasFWHM()
		#print(' total area after = ',self.total_surface_area)
		self.plotting_table()
	
	
	def SaveFigureAll(self): # save figure any format 
		self.file_to_save_fig = asksaveasfilename(initialdir="", \
		filetypes =(("PDF File","*.pdf"),("PNG File", "*.png"),("SVG File","*.svg")),
		   title = "File to name.")
		print(' file figure ',self.file_to_save_fig)
		self.fig.savefig(self.file_to_save_fig, dip=600,bbox_inches='tight')
		
		
		
		
	
	
	
	
	
	def OpenParameterFile(self):
		self.input_fit_file = askopenfilename(initialdir="",
								   filetypes =(("All Files","*.*"),("Text File", "*.txt"),("Data File","*.dat"),("Input File","*.inp")),
								   title = "Choose a file."
								   )
		#print (self.input_fit_file)
		#Using try in case user types in unknown file or closes without choosing a file.
		try:
			self.User_ParmGuess = []
			self.User_Parm_Flag = []
			self.User_FuncType  = []
			with open(self.input_fit_file,'r') as FitGuessFile:
				for line in FitGuessFile.readlines():
					data = line.split()
					#print( data)
					if( data[0][0] != '#'):
						self.User_FuncType.append (        data[0]  )  # type of function 'Hubbert' 'Gauss' or 'Lorentz'
						self.User_ParmGuess.append( float( data[2] ))  # parameter A 
						self.User_ParmGuess.append( float( data[3] ))  # parameter B 
						self.User_ParmGuess.append( float( data[4] ))  # parameter C
						function_flag_read =( 1  if data[5]=='fix' else 0)
						self.User_Parm_Flag.append( function_flag_read  )  
						# fix or opt   fix - parameters are kept fixed, opt - they are optimized
			
			#print('rozmiar ', len(self.User_FuncType))
			#print('print Parm Flag  and Func Type ', self.User_Parm_Flag, self.User_FuncType)			
			#print(' input par.')
			for i in range(len(self.User_FuncType)):
				#print(' {} {} {} {} {}'.format(self.User_FuncType[i],self.User_ParmGuess[i],self.User_ParmGuess[i+1],\
				#self.User_ParmGuess[i+2],self.User_Parm_Flag[i]))
				self.number_of_peaks = len(self.User_FuncType)
				
			#self.Show_Parameters_GUI()	
			if  len(self.User_FuncType)==0 :
				raise ValueError
				
		except IOError:
			messagebox.showinfo("ERROR", "No file exists or is empty")
		except ValueError :
			messagebox.showerror("ERROR", "something wrong with the size ")
		else:
			if    len(self.WidgetFixOptChecks) == 0:
				self.CreateParameterTable()
			elif  len(self.WidgetFixOptChecks) > 0:	
				if self.number_of_parameters == self.number_of_peaks_in_table:
					pass #update table only
				else:	
					self.RemoveParameterTable()
					self.CreateParameterTable()
			self.UpdateParameterTable(self.User_FuncType, self.User_Parm_Flag, self.User_ParmGuess)
			self.labelsBottom.configure(text='Component functions: read from file')
			self.plotting_table()
			#self.plotting_updating()		
			#print("No file exists")

	def UpdateParameterTable(self, Function_Type, Optimization_Flag, Parameter_Value  ):
		# don't care about table sizes and number of parameters, table asssumed to exist and size = number of peaks  
		#print('\n\n\n\niside updating parameter table')
		#print('types ',Function_Type)
		#print('flags ',Optimization_Flag)
		#print('param ',Parameter_Value)
		#print('\n'*4)
		for i in range(self.number_of_peaks):
			self.WidgetFunType_Var[i].set( Function_Type[i] )
			self.WidgetFixOptChecks_Var[i].set( Optimization_Flag[i] )
			
			#print('flag ',self.WidgetFixOptChecks_Var[i].get())
			for j in range(self.number_of_parameters):
				self.WidgetParameters_Var[i*self.number_of_parameters +j].set( Parameter_Value[i*self.number_of_parameters +j])
				self.WidgetParameters_Obj[i*self.number_of_parameters +j].configure(format="%0.5f")
				
			
			
				
	
	def OpenDataFile(self):
		
		self.input_data_file = askopenfilename(initialdir="",
								   filetypes =(("All Files","*.*"),("Text File", "*.txt"),("Data File","*.dat"),("Fit Files","*.fit")),
								   title = "Choose a file."
								   )
								
								
		#print (self.input_data_file)
		#Using try in case user types in unknown file or closes without choosing a file.
		if( len(self.input_data_file)>0):
			self.xdata =[]
			self.ydata =[]
			self.ydata_org =[]
			try:	
				with open (self.input_data_file,'r') as f:
					for line in f.readlines():
						data = line.split()
						self.xdata.append( float( data[0] ))
						self.ydata.append( float( data[1] ))
						self.ydata_org.append( float(data [1]))
					self.MinMaxDataValues =[ min(self.xdata), max(self.xdata), min(self.ydata_org), max(self.ydata_org)]
				
			except FileNotFoundError:
				error_message='File with data to fit ({}) does not exist'.format(name); 
				messagebox.showinfo("ERROR", error_message)
			self.data_range  =[min(self.xdata),max(self.xdata),min(self.ydata_org),max(self.ydata_org)]
			#print(' data range =',self.data_range)
			self.plot_range  =[]
			for j in range(4):	self.plot_range.append(self.data_range[j]) 
			#print('data range = ',self.data_range )
			#print('plot range = ',self.plot_range )
			self.AutoscalePlotRange()
			self.plotting_data()					
			
			
			self.x_at_maximal_y_value = 0
			for j in range(len(self.xdata)):
				if self.ydata_org[j] == self.data_range[3]: self.x_at_maximal_y_value = self.xdata[j]
			#print(' x at maximal y is ',self.x_at_maximal_y_value)
			
			
			
			#self.InitialGuess 		= [] # initial guess if none is provided by User 
			self.DistributionRange = []  # range of x and y values that contains majority of data 
			self.GuessScale        = []  # scale step for n-initial functions to guess
			
			
			self.DistributionRange =[max(self.xdata),0,min(self.xdata),0] # x,y (left) and xy (right) pair of values x,y for y=f(x) 
			for i in range(len(self.ydata)):
				#if ( self.ydata[i]>0.15*max(self.ydata)  and self.ydata[i]<0.25*max(self.ydata) ): 
				if ( self.ydata[i]>0.05*max(self.ydata)  and self.ydata[i]<0.10*max(self.ydata) ): 
					# search among the y=f(x) for which y is in (0.15;0.25) of y_max
					if  self.xdata[i] < self.DistributionRange[0] :
						self.DistributionRange[0]=self.xdata[i]
						self.DistributionRange[1]=self.ydata[i]
					elif self.xdata[i] > self.DistributionRange[2]:
						self.DistributionRange[2]=self.xdata[i]
						self.DistributionRange[3]=self.ydata[i]
			#print(' distribution range (x,y (left) and x,y (right)',self.DistributionRange)
			if  len( self.xdata ) > 0:
				xmax_spin = self.data_range[1]
				ymax_spin = self.data_range[3]
				wmax_spin = abs(self.data_range[1] - self.data_range[0])
			else:
				ymax_spin = 1.0
				xmax_spin = 100.0
					
			for i in range( self.number_of_peaks):
				self.WidgetParameters_Obj[1+i*3].configure(to=xmax_spin)
				self.WidgetParameters_Obj[  i*3].configure(to=ymax_spin)
				self.WidgetParameters_Obj[2+i*3].configure(to=wmax_spin)
				
					

				
			self.GenerateInitialGuess()
			self.UpdateParameterTable( ['Gauss']*self.number_of_peaks , [0]*self.number_of_peaks, self.InitialGuess)
			self.labelsBottom.configure(text='Component functions: initial guess')
						
			
	
	def calculate_function(self):
		if( self.selected_function_type.get() == 'Gauss'):
				self.outputBox.delete(0,len(self.result1.get()))
				self.outputBox.configure(foreground='red')
				self.outputBox.insert(0,'{}'.format(12.2))
				self.outputBox.configure(state='disabled')
				
		
	def _quit(self):
		self.master.quit()
		self.master.destroy()
		exit() 	
	
	
	def ParametersChangedinTable(self):
		# read only affected 
		pass 
		
	
	
	def RemoveParameterTable(self):
		# distroy widgets 
		#print(' length ',len(self.WidgetFixOptChecks))
		if( len(self.WidgetFixOptChecks) >0): # check if the table was created before 
			for child in self.labelsBottom.winfo_children():
				child.destroy()
		# clear all table lists 
		self.WidgetParameters_Obj   =[]
		self.WidgetFixOptChecks     =[]
		self.WidgetFixOptChecks_Var =[]
		self.WidgetPlotChecks_Var   =[]
		self.WidgetParameters_Var   =[]
		self.WidgetPlotChecks       =[]
		self.WidgetLabel		     =[]
		self.WidgetFunType_Obj      =[]
		self.WidgetFunType_Var      =[]
		self.WidgetPlotRemChecks_Var=[]
		self.WidgetIncremenet_Var =[]
		
				
		
	def UpdateGlobalPlotingComponents(self):
		if self.CheckShowComponentsFit_Var.get() == False:
			for i in range(self.number_of_peaks):
				self.WidgetPlotChecks_Var[i].set(False)
		else:
			for i in range(self.number_of_peaks):
				self.WidgetPlotChecks_Var[i].set(True)
		self.plotting_table()		
		
			
		
	
	
	def CreateParameterTable(self): 
		'''self.AddPeakButton =  ttk.Button(self.labelsBottom, text='+',width=2)
		self.AddPeakButton.grid( column=9,row=6)
		self.RemovePeakButton= ttk.Button(self.labelsBottom, text='-',width=2)
		self.RemovePeakButton.grid(column=10,row=6)'''
		self.RemoveParameterTable()
		ttk.Label(self.labelsBottom, text='#',foreground='tomato',justify=tk.RIGHT).grid(column=0,row=0,padx=2)
		ttk.Label(self.labelsBottom, text='Type',foreground='tomato').grid(column=1,row=0)
		ttk.Label(self.labelsBottom, text='Parameter (1)',foreground='tomato').grid(column=2,row=0)
		ttk.Label(self.labelsBottom, text='Parameter (2)',foreground='tomato').grid(column=3,row=0)
		ttk.Label(self.labelsBottom, text='Parameter (3)',foreground='tomato').grid(column=4,row=0)
		
		ttk.Label(self.labelsBottom, text='Fix',foreground='tomato',justify=tk.LEFT).grid(column=5,row=0)
		ttk.Label(self.labelsBottom, text='Plot',foreground='tomato').grid(column=6,row=0)
		ttk.Label(self.labelsBottom, text='Del',foreground='tomato').grid(column=7,row=0)
		ttk.Label(self.labelsBottom, text='Area',foreground='tomato').grid(column=8,row=0)
		ttk.Label(self.labelsBottom, text='%',foreground='tomato').grid(column=9,row=0)
		ttk.Label(self.labelsBottom, text='FWHM',foreground='tomato').grid(column=10,row=0,padx=1)
		
		#print('number of peaks = ',self.number_of_peaks )
		self.WidgetParameters_Obj   	=[]
		self.WidgetFixOptChecks    		=[]
		self.WidgetFixOptChecks_Var 	=[]
		self.WidgetPlotChecks_Var   	=[]
		self.WidgetParameters_Var   	=[]
		self.WidgetPlotChecks       	=[]
		self.WidgetLabel		     	=[]
		self.WidgetFunType_Obj      	=[]
		self.WidgetFunType_Var      	=[]
		self.WidgetPlotChecks			=[]
		self.WidgetPlotChecks_Var		=[]
		self.WidgetIncremenet_Var		=[]
		self.WidgetIncrementEntry		=[]
		self.WidgetComponentArea        =[] 
		self.WidgetComponentWeight		=[]
		
		
		for i in range( self.number_of_peaks):
			self.WidgetFixOptChecks.append('')	# CheckBox for Fit/Only Plot flag (0/1)
			self.WidgetFixOptChecks_Var.append ( tk.IntVar()		)	# 
			self.WidgetPlotChecks.append('')
			self.WidgetPlotChecks_Var.append( tk.BooleanVar())
			self.WidgetLabel.append('')
			self.WidgetPlotChecks_Var.append( tk.BooleanVar())
			self.WidgetFunType_Obj.append('')
			self.WidgetFunType_Var.append( tk.StringVar())
			self.WidgetComponentArea.append('')
			self.WidgetComponentWeight.append('')
			self.WidgetComponentFWHM.append('')
			self.WidgetPlotRemChecks_Var.append( tk.IntVar() )
			self.WidgetPlotRemChecks.append('')
		
		for i in range(self.number_of_parameters):	
			self.WidgetIncrementEntry.append('WP_Increment')
			self.WidgetIncremenet_Var.append( tk.StringVar() )
			
		for i in range( self.number_of_parameters * self.number_of_peaks):	
			self.WidgetParameters_Obj.append('WP_'		+str(i)	)
			self.WidgetParameters_Var.append(tk.DoubleVar())	

		for i in range( self.number_of_peaks):
			for j in range(self.number_of_parameters):
				self.WidgetParameters_Obj[j+i*3]=tk.Spinbox(self.labelsBottom,width=10,from_=0.0001,to=100.0,\
				increment=self.ParameterIncrease[j],format="%10.5f",textvariable=self.WidgetParameters_Var[j+i*3],\
				justify=tk.RIGHT,command=self.plotting_table_and_update)
				
				self.WidgetParameters_Obj[3*i+j].grid(column=2+j, row=1+i)
				self.WidgetParameters_Obj[j+i*3].bind("<Return>",self.UpdateParametersUser)
				#self.WidgetParameters_Var[j+i*3].trace('w',self.UpdateParametersUser)
		
		if  len( self.xdata ) > 0:
			xmax_spin = self.data_range[1]
			ymax_spin = self.data_range[3]
			wmax_spin = abs(self.data_range[1] - self.data_range[0])
		else:
			ymax_spin = 1.0
			xmax_spin = 100.0
			wmax_spin = 100.0
				
		for i in range( self.number_of_peaks):
			self.WidgetParameters_Obj[1+i*3].configure(to=xmax_spin)
			self.WidgetParameters_Obj[  i*3].configure(to=ymax_spin)
			self.WidgetParameters_Obj[2+i*3].configure(to=wmax_spin)
				
				
			
		for i in range(self.number_of_peaks): # fix opt flag checkbox
			self.WidgetComponentArea[i]=ttk.Label(self.labelsBottom,text='surf.area',foreground='sea green')
			self.WidgetComponentArea[i].grid(column=8,row=i+1,padx=5)
			self.WidgetComponentWeight [i]=ttk.Label(self.labelsBottom,text='weight%',foreground='salmon')
			self.WidgetComponentWeight[i].grid(column=9,row=i+1,padx=5)
			self.WidgetComponentFWHM[i]=ttk.Label(self.labelsBottom,text='FWHM',foreground='salmon')
			self.WidgetComponentFWHM[i].grid(column=10,row=i+1,padx=5)
	
			# number of function 
			self.WidgetLabel[i]=ttk.Label(self.labelsBottom,text='{:2d}'.format(i+1),foreground='salmon')
			self.WidgetLabel[i].grid(column=0,row=1+i,sticky=tk.W,padx=2,pady=1)
			#self.WidgetLabel[i].configure(foreground= self.PlotsColorList[i])
			
			# checkBox for optimize or fixed parameter value 
			self.WidgetFixOptChecks[i]=tk.Checkbutton(self.labelsBottom,text='',\
			variable=self.WidgetFixOptChecks_Var[i],command=self.check_opt_fix_flags)
			self.WidgetFixOptChecks[i].grid(column=5,row=1+i,sticky=tk.W,pady=1)
			self.WidgetFixOptChecks[i].deselect()
			
			
			# plot (show) or not 
			self.WidgetPlotChecks[i]=tk.Checkbutton(self.labelsBottom,text='',variable=self.WidgetPlotChecks_Var[i],\
			command= lambda i=i: self.ShowHide(i))
			#self.WidgetPlotChecks[i].grid(column=7,row=6+i,sticky=tk.W,padx=1,pady=1)
			self.WidgetPlotChecks[i].grid(column=6,row=1+i,sticky=tk.W,pady=1)
			self.WidgetPlotChecks[i].select()
			
			# function type 
			self.WidgetFunType_Obj[i]=ttk.Combobox(self.labelsBottom,width=8,textvariable=self.WidgetFunType_Var[i],\
			state='readonly',justify=tk.CENTER)
			self.WidgetFunType_Obj[i]['values']=('Gauss','Hubbert','Lorentz')
			self.WidgetFunType_Obj[i].current(0)
			self.WidgetFunType_Obj[i].grid(column=1,row=1+i,sticky=tk.W,padx=2,pady=2)
			
			# remove this function 
			self.WidgetPlotRemChecks[i]=tk.Checkbutton(self.labelsBottom,text='')#,variable=self.WidgetPlotRemChecks_Var[i])
			self.WidgetPlotRemChecks[i].grid(column=7,row=1+i,sticky=tk.W,pady=1)
			self.WidgetPlotRemChecks[i].deselect()
		
		# increment in table (last row)
		ttk.Label(self.labelsBottom,text='Increment =',justify=tk.RIGHT,foreground='steel blue').grid(column=1,row=self.number_of_peaks+1)
		for i in range(self.number_of_parameters):
			self.WidgetIncrementEntry[i]=tk.Entry(self.labelsBottom, textvariable=self.WidgetIncremenet_Var[i],foreground='steel blue', \
			width=6, justify=tk.RIGHT)
			self.WidgetIncremenet_Var[i].set("0.01")
			self.WidgetIncrementEntry[i].grid(column=2+i,row=self.number_of_peaks+1,sticky=tk.E, padx=5,pady=5)
			self.WidgetIncremenet_Var[i].trace('w',self.UpdateIncrement)
		
		
		self.CheckShowComponentsFit_Var = tk.BooleanVar()
		self.CheckShowComponentsFit_Obj = tk.Checkbutton(self.labelsBottom,variable=self.CheckShowComponentsFit_Var,text='plot components',\
		command=self.UpdateGlobalPlotingComponents)
		self.CheckShowComponentsFit_Obj.grid(column=6,row=self.number_of_peaks+1,columnspan=4,sticky=tk.W,padx=1,pady=2)
		self.CheckShowComponentsFit_Obj.select()
		
		self.GenerateInitialGuess()
		# here generate accordingly 
		self.InitialTypes = [self.selected_function_type.get()]*self.number_of_peaks
		self.InitialFlags = ['0']*self.number_of_peaks
		
		# here update with the parameter values if any existed 
		if len(self.Current_Function_Types) > 0:
			if len(self.Current_Function_Types) >= self.number_of_peaks:
				for i in range(self.number_of_peaks):
					self.InitialFlags[i]= self.Current_Function_Flags[i]
					self.InitialTypes[i]= self.Current_Function_Types[i]
					for j in range(self.number_of_parameters):
						self.InitialGuess[i*3 +j] = self.Current_Parameters[i*3 +j]
			elif len(self.Current_Function_Types) < self.number_of_peaks:
				for i in range(len(self.Current_Function_Types)):
					self.InitialFlags[i]= self.Current_Function_Flags[i]
					self.InitialTypes[i]= self.Current_Function_Types[i]
					for j in range(self.number_of_parameters):
						self.InitialGuess[i*3 +j] = self.Current_Parameters[i*3 +j]
		self.UpdateParameterTable( self.InitialTypes, self.InitialFlags, self.InitialGuess)	
				
			
	
	
	
	
	def check_opt_fix_flags(self):
		for i in range(self.number_of_peaks):
			if( self.WidgetFixOptChecks_Var[i].get() == 1):
				self.WidgetFunType_Obj[i].configure(state='disable')
				self.WidgetLabel[i].configure(foreground='grey')
				for j in range(3):
					self.WidgetParameters_Obj[i*self.number_of_parameters + j].configure(state='disable')
			if( self.WidgetFixOptChecks_Var[i].get() == 0):
				self.WidgetFunType_Obj[i].configure(state='normal')
				self.WidgetLabel[i].configure(foreground='tomato')
				for j in range(3):
					self.WidgetParameters_Obj[i*self.number_of_parameters + j].configure(state='normal')		
				
				
	
	
	
			
	def plotting_data(self):
		if len(self.xdata) > 0:
			if self.CheckShowDataPlot_Var.get() == True: 
				self.line1.set_ydata([xx for xx in self.ydata_org])
				self.line1.set_xdata([xx for xx in self.xdata])
				self.ax.axis([self.plot_range[0],self.plot_range[1],self.plot_range[2],self.plot_range[3]])
				self.canvas.draw()
			else:
				self.line1.set_ydata([])
				self.line1.set_xdata([])
				self.ax.axis([self.plot_range[0],self.plot_range[1],self.plot_range[2],self.plot_range[3]])
				self.canvas.draw()
			
				
	
			
				
	


# ====================================================================
#  MAIN PROGRAM LOOP
# ====================================================================
random.seed()
root = tk.Tk()
root.update()
root.minsize(root.winfo_width(), root.winfo_height())
app = Application(master=root)
root.update()
root.minsize(root.winfo_width(), root.winfo_height())
app.mainloop()
	
	
	