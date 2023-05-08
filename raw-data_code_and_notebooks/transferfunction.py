#!/usr/bin/env python


import numpy as np
from scipy.optimize import curve_fit, minimize, leastsq


class TransferFunction():

    def __init__(self, inputs, outputs, name):

        self.name = name
        self.inputs = inputs #system inputs, in RPUs preferably
        self.outputs = outputs #measured outputs
        self.ymax = max(outputs)
        # self.ymin = None # AJT added
        self.ymin = min(outputs)/3 # AJT commented out
        self.kd = None
        self.n = None

    def fit_data(self, cost_function='leastsq', x_initial=[0.5,2,0.001]): # AJT added third value of 0.001
        #fit Kd and n, use initial guesses of 0.5 and 2, which are pretty good
        #AJT: had been using 10, 2 for activator function and looks fine for repressor function except when input is RPU
        #AJT: Pbad requires initial guess of higher (100,2) because Kd is so high (requires a lot of inducer)
        
        # set bounds for x_initial values ((Kdlow,Kdhigh),(nlow,nhigh),(yminlow,yminhigh))
        # bnds = ((0,1),(0,20),(min(0.001,min(self.outputs)/10),min(self.outputs))) # AJT added this line or next
        # bnds = ((0,1),(0,20),(min(self.outputs)/3,min(self.outputs)))

        res = minimize(self._cost, x0=x_initial, args=(cost_function), method='Nelder-Mead') # AJT commented out
        # res = minimize(self._cost, x0=x_initial, args=(cost_function), method='TNC', bounds=bnds) # AJT added since 'Nelder-Mead' method doesn't allow bounds
        
        #elif cost_function == 'l1':
        self.kd = res.x[0]
        self.n = res.x[1]
        # self.ymin = res.x[2] # AJT added
        return

    def _cost(self, parameters, cost_function): 
    
        #minimize function can take only 1 parameter, can pass other as args after later on
        #should take as parameters K and n, which we vary
        #plug them into the hill equation, calculate some error
        #then minimize the error

        #x is a array of input points, y is the measured output points
        #then compute the ypred and error based on parameters
        k = parameters[0]
        n = parameters[1]
        # ymin = parameters[2] # AJT added
        x = np.array(self.inputs)
        y = np.array(self.outputs)
        
        ypred = self._hill_eqn(x, self.ymax, self.ymin, k, n) # AJT changed self.ymin to ymin
        #ypred = hill_eqn1(x, ymax, ymin, k, n)
        #ypred = hill_eqn(x, ymax, ymin, k, n)
        
        if cost_function == 'leastsq':
        	err = sum( ((y - ypred)/y)**2) #standard least squares but normalized
        elif cost_function == 'l1':
        	err = sum( (np.absolute(y - ypred))/y) #l1 normalized
        else:
        	raise TypeError("Unknown cost function. Options are 'leastsq' or 'l1'.")
        
        return err
    
    def _hill_eqn(self, x, ymax, ymin, K, n):
    	#internal hill eqn for used for fitting data
        return ymin + (ymax - ymin)*(K**n)/( (x**n) + (K**n))
        
    def hill_eqn(self, x):
    	'''Given 'x', calculates output via the following hill equation:
    	y = ymin + (ymax - ymin) * kd^n / (x^n + kd^n). Returns a single y value'''
    	
    	assert self.kd is not None, 'TransferFunction kd attribute is not set'
    	assert self.n is not None, 'TransferFunction n attribute is not set'
    	assert self.ymin is not None, 'TransferFunction ymin attribute is not set' # AJT added
    	return self.ymin + (self.ymax - self.ymin)*(self.kd**self.n)/( (x**self.n) + (self.kd**self.n)) # k**self.n should be numerator for repression function

def build_transfer_functions(inputs, outputs, names):
	'''Convenience function to bulk generate a bunch of TransferFunction objects.
	'inputs' = list of inputs values (must be the same for each TransferFunction)
	'outputs' = list of output value lists = [[y1,y2,...], ...]. sublists must be 
	same size as 'inputs'. 
	'names' = list of names. 
	Returns a list of TransferFunction objects.'''
	
	
	# for output in outputs:
	# 	if len(inputs) != len(output):
	# 		raise TypeError('Input and output arrays have different lengths.')
	
	if len(outputs) != len(names):
		raise TypeError('Must specify same number of names as output sublists')
		
	transfer_functions = []
	for i in range(len(outputs)):
		tf = TransferFunction(inputs[i], outputs[i], names[i])
		tf.fit_data()
		transfer_functions.append(tf)
		
	return transfer_functions
	    
	
if __name__ == "__main__":
	main()	
	

	
	
	