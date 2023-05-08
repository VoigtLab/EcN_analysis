#!/usr/bin/env python


import numpy as np
from scipy.optimize import curve_fit, minimize, leastsq
import warnings


class TransferFunction():

    def __init__(self, inputs, outputs, name, kind):#, x_units):

        self.name = name
        self.inputs = inputs #system inputs, in RPUs preferably
        self.outputs = outputs #measured outputs
        self.kind = kind
        # self.ymax = max(outputs)
        # self.ymin = None # AJT added
        if self.kind == 'activator':
            self.ymin = min(outputs) # AJT commented out and added /2
            self.ymax = None # AJT added
            ## AJT added below
            # x_init = [kd, n, ymin]
            if max(outputs) > 10: # for when y-axis is in au
                self.x_init = [0.5,2,1000]
            else:
                self.x_init = [0.5,2,1]
        elif self.kind == 'repressor': # AJT added
            # self.ymin = min(outputs) # AJT added
            self.ymin = None
            self.ymax = max(outputs) # AJT added
            # AJT added below
            # x_init = [kd, n, ymin]
            if max(outputs) > 10: # for when y-axis is in au
                self.x_init = [0.5,2,50]
            else:
                self.x_init = [0.5,2,0.01]
        # self.x_units = x_units
        self.kd = None
        self.n = None

    # def fit_data(self, cost_function='leastsq', x_initial=[0.5,2]): # AJT added third value of 0.001
    # def x_initial(self)
    #     # # AJT added
    #     if self.kind == 'activator':
            

    #     if self.kind == 'repressor':
    #         if max(outputs) > 10: # for when y-axis is in au
    #             x_init = [0.5,2,50]
    #         else:
    #             x_init = [0.5,2,0.01]


    def fit_data(self, cost_function='leastsq'):#, x_initial=self.x_init: # for allowing ymax to float, AJT added third value of 
    # def fit_data(self, cost_function='leastsq', x_initial=[0.5,2,0.001]): # for allowing ymax to float, AJT added third value of 
        #fit Kd and n, use initial guesses of 0.5 and 2, which are pretty good
        #AJT: had been using 10, 2 for activator function and looks fine for repressor function except when input is RPU
        #AJT: Pbad requires initial guess of higher (100,2) because Kd is so high (requires a lot of inducer)
        
        # set bounds for x_initial values ((Kdlow,Kdhigh),(nlow,nhigh),(yminlow,yminhigh))
        # bnds = ((0,1),(0,20),(min(0.001,min(self.outputs)/10),min(self.outputs)))
        # bnds = ((0,1),(0,20),(min(self.outputs)/3,min(self.outputs)))
        bnds_rep = ((min(self.inputs),max(self.inputs)),(0,20),(min(self.outputs)/3,min(self.outputs)))

        # set bounds for x_initial values ((Kdlow,Kdhigh),(nlow,nhigh),(ymaxlow,ymaxhigh))
        # bnds = ((0,1),(0,20),(max(0.001,max(self.outputs)/10),max(self.outputs)))
        # bnds_act = ((0.5,0.9),(2.5,3),(max(self.outputs)*1.25,max(self.outputs)*2))
        bnds_act = ((min(self.inputs),max(self.inputs)),(0,20),(max(self.outputs),max(self.outputs)))
        
        

        if self.kind == 'repressor':
            # res = minimize(self._cost, x0=x_initial, args=(cost_function), method='Nelder-Mead') # AJT commented out
            # res = minimize(self._cost, x0=x_initial, args=(cost_function), method='TNC', bounds=bnds_rep) # AJT added since 'Nelder-Mead' method doesn't allow bounds
            res = minimize(self._cost, x0=self.x_init, args=(cost_function), method='TNC', bounds=bnds_rep) # AJT added since 'Nelder-Mead' method doesn't allow bounds

        if self.kind == 'activator':
            # res = minimize(self._cost, x0=x_initial, args=(cost_function), method='TNC', bounds=bnds_act) # AJT added since 'Nelder-Mead' method doesn't allow bounds
            res = minimize(self._cost, x0=self.x_init, args=(cost_function), method='TNC', bounds=bnds_act) # AJT added since 'Nelder-Mead' method doesn't allow bounds
        
        #elif cost_function == 'l1':
        self.kd = res.x[0]
        self.n = res.x[1]
        # self.ymin = res.x[2] # AJT added
        if self.kind == 'activator':
            self.ymax = res.x[2] # AJT added

        if self.kind == 'repressor':
            self.ymin = res.x[2] # AJT added

        # cof = np.reshape(np.array(res.x),(-1,1))
        # print(res.x, cof)
        # self.r_squared = 1 - np.square(np.array(self.outputs)-np.array(self.inputs).dot(cof)).sum() / (np.var(np.array(self.outputs)) * len(self.outputs))

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
        if self.kind == 'activator':
            ymax = parameters[2] # AJT added

        if self.kind == 'repressor':
            ymin = parameters[2] # AJT added

        x = np.array(self.inputs)
        y = np.array(self.outputs)
        
        if self.kind == 'repressor':
            ypred = self._hill_eqn_rep(x, self.ymax, ymin, k, n) # AJT changed self.ymin to ymin
        elif self.kind == 'activator':
            # ypred = self._hill_eqn_act(x, self.ymax, self.ymin, k, n) # AJT changed self.ymin to ymin
            ypred = self._hill_eqn_act(x, ymax, self.ymin, k, n) # AJT changed self.ymax to ymax
        else:
            warnings.warn("_cost does not recognize kind")
        #ypred = hill_eqn1(x, ymax, ymin, k, n)
        #ypred = hill_eqn(x, ymax, ymin, k, n)
        
        if cost_function == 'leastsq':
            err = sum( ((y - ypred)/y)**2) #standard least squares but normalized
        elif cost_function == 'l1':
            err = sum( (np.absolute(y - ypred))/y) #l1 normalized
        else:
            raise TypeError("Unknown cost function. Options are 'leastsq' or 'l1'.")
        
        return err

    def _hill_eqn_rep(self, x, ymax, ymin, K, n):
        #internal hill eqn for used for fitting data
        return ymin + (ymax - ymin)*(K**n)/( (x**n) + (K**n))
        
    def hill_eqn_rep(self, x):
        '''Given 'x', calculates output via the following hill equation:
        y = ymin + (ymax - ymin) * kd^n / (x^n + kd^n). Returns a single y value'''
        
        assert self.kd is not None, 'TransferFunction kd attribute is not set'
        assert self.n is not None, 'TransferFunction n attribute is not set'
        assert self.ymin is not None, 'TransferFunction ymin attribute is not set' # AJT added
        return self.ymin + (self.ymax - self.ymin)*(self.kd**self.n)/( (x**self.n) + (self.kd**self.n)) # k**self.n should be numerator for repression function

    def _hill_eqn_act(self, x, ymax, ymin, K, n):
        #internal hill eqn for used for fitting data
        return ymin + (ymax - ymin)*(x**n)/( (x**n) + (K**n))
        
    def hill_eqn_act(self, x):
        '''Given 'x', calculates output via the following hill equation:
        y = ymin + (ymax - ymin) * kd^n / (x^n + kd^n). Returns a single y value'''
        
        assert self.kd is not None, 'TransferFunction kd attribute is not set'
        assert self.n is not None, 'TransferFunction n attribute is not set'
        # assert self.ymin is not None, 'TransferFunction ymin attribute is not set' # AJT added
        assert self.ymax is not None, 'TransferFunction ymax attribute is not set' # AJT added
        return self.ymin + (self.ymax - self.ymin)*(x**self.n)/( (x**self.n) + (self.kd**self.n)) # k**self.n should be numerator for repression function

def build_transfer_functions(inputs, outputs, names, kind):
    '''Convenience function to bulk generate a bunch of TransferFunction objects.
    'inputs' = list of inputs values (must be the same for each TransferFunction)
    'outputs' = list of output value lists = [[y1,y2,...], ...]. sublists must be 
    same size as 'inputs'. 
    'names' = list of names. 
    Returns a list of TransferFunction objects.'''


    # for output in outputs:
    #   if len(inputs) != len(output):
    #       raise TypeError('Input and output arrays have different lengths.')

    if len(outputs) != len(names):
        raise TypeError('Must specify same number of names as output sublists')

    if kind not in ['repressor','activator']:
        raise TypeError("Unknown kind. Options are 'activator' or 'repressor'")

    transfer_functions = []
    for i in range(len(outputs)):
        tf = TransferFunction(inputs[i], outputs[i], names[i], kind)
        tf.fit_data()
        transfer_functions.append(tf)
        
    return transfer_functions


if __name__ == "__main__":
    main()  



