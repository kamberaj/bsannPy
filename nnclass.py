import numpy as np
import sys
#---------------------------------------------------------------------------------
#### Miscellaneous functions
def square(x):
    return x*x

def TransferFunctionInn(z, ifunc):
    """The transfer function."""
    if ifunc == 1:
       return 1.0/(1.0 + np.exp(-z))
    elif ifunc == 2:
       return 1.0/(1.0 + square(z))
    else:
       return np.tanh(z)

def TransferFunctionInn_derivative(z, ifunc):
    """Derivative of the transfer function."""
    if ifunc == 1:
       return TransferFunctionInn(z,ifunc) * (1.0 - TransferFunctionInn(z,ifunc))
    elif ifunc == 2:
       return ( -2.0 * z * square(TransferFunctionInn(z,ifunc)) )
    else:
       return 1.0 - square( TransferFunctionInn(z,ifunc) )

#-----------------------------------------------------------------------------------
class Config:
    # Gradient descent parameters
    Weight = 0.0001        # learning rate for gradient descent
    C1     = 0.01          # regularization strength
    C2     = 0.01          # Individual knowledge strength
    C3     = 0.01          # Individual knowledge strength
#--------------------------------------------------------------------------------
class Layer:
    def __init__(self,data_size,dims,nneurons, w, b):
        self.nneurons = nneurons
        self.w = w
        self.b = b
        self.w_opt = w
        self.b_opt = b
        self.w_global_opt = w
        self.b_global_opt = b
        self.dw = np.zeros(w.shape)
        self.db = np.zeros(b.shape)
        self.dims = dims 
        self.data_size = data_size
        self.output = np.zeros((data_size,dims), float)

    def setCurrentParams(self, w, b):
        self.w = w
        self.b = b

    def getCurrentParams(self):
        return self.w, self.b

    def setOptimalParams(self, w, b):
        self.w_opt = w
        self.b_opt = b

    def getOptimalParams(self):
        return self.w_opt, self.b_opt

    def setGlobalBestParams(self, w, b):
        self.w_global_opt = w
        self.b_global_opt = b

    def getGlobalBestParams(self):
        return self.w_global_opt, self.b_global_opt

#--------------------------------------------------------------------------------
class NeuralNetwork:
    def __init__(self, x, y, dims, LayerNeurons, ifunc):
        self.input      = x
        self.y          = y
        self.data_size  = x.shape[0]
        self.nlayers    = len(LayerNeurons)-1
        self.output     = np.zeros(self.y.shape)
        self.output_opt = np.zeros(self.y.shape)
        self.err_best   = 1000000.0
        self.layer = []
        self.TransferFunction = ifunc

        for i in range(self.nlayers):
            w = np.random.randn(LayerNeurons[i], LayerNeurons[i+1]) / np.sqrt(LayerNeurons[i])
            b = np.zeros( LayerNeurons[i+1], float)
            self.layer.append(Layer(self.input.shape[0], dims[i], LayerNeurons[i], w, b))

##########################################################################################################    
    def NNrun(self, epochs=100): 
        for i in range(epochs):
            self.feedforward()
            err, best = self.getOptimalNetwork()
            self.backprop()

        self.feedforward()
        err, best = self.getOptimalNetwork()
##########################################################################################################
    def getOptimalNetwork(self):
        err = self.rmse(self.output)
        if err < self.err_best:
           self.err_best = err
           self.output_opt = self.output
           for k in range(self.nlayers):
               self.layer[k].setOptimalParams( self.layer[k].w, self.layer[k].b )
        return err, self.err_best
##########################################################################################################
    def predict(self, x, myLayer):
        self.activation  = x
        self.activations = [x]
        self.zs = []
        for i in range(self.nlayers-1):
            b = np.zeros( (x.shape[0],self.layer[i+1].nneurons), float)
            for i1 in range(x.shape[0]):
                for j1 in range(self.layer[i+1].nneurons):
                    b[i1][j1] = self.layer[i].b_opt[j1]
            
            z = np.dot(self.activation,  self.layer[i].w_opt) + b
            self.zs.append(z)
            self.activation = TransferFunctionInn(z,self.TransferFunction)
            self.activations.append(self.activation)
       
        z = np.dot(self.activations[-1], self.layer[self.nlayers-1].w_opt) + self.layer[self.nlayers-1].b_opt
        self.zs.append(z)
        self.activation = TransferFunctionInn(z,self.TransferFunction)
        self.activations.append(self.activation)
 
        return self.activations[-myLayer]

##########################################################################################################
    def getActivation(self, k):
        return self.activations[-k]
##########################################################################################################            
    def feedforward(self):  #forward feed
        self.activation  = self.input
        self.activations = [self.input]
        self.zs = []
        for i in range(self.nlayers-1):
            b = np.zeros( (self.activation.shape[0], self.layer[i+1].nneurons), float)
            for i1 in range(self.activation.shape[0]):
                for j1 in range(self.layer[i+1].nneurons):
                    b[i1][j1] = self.layer[i].b[j1]
           
            z = np.dot(self.activation,  self.layer[i].w) + b
            self.zs.append(z)
            self.activation = TransferFunctionInn(z,self.TransferFunction)
            self.activations.append(self.activation)
        
        z = np.dot(self.activations[-1], self.layer[self.nlayers-1].w) + self.layer[self.nlayers-1].b
        self.zs.append(z)
        self.activation = TransferFunctionInn(z,self.TransferFunction)
        self.activations.append(self.activation)
        self.output = self.activations[-1]
##########################################################################################################
    def backprop(self):   # backward pass
        delta = 2.0*self.cost_derivative(self.activations[-1], self.y) * TransferFunctionInn_derivative(self.zs[-1],self.TransferFunction)
        self.layer[-1].db = np.sum( delta, 0 )
        self.layer[-1].dw = np.dot(self.activations[-2].T, delta)
        for k in xrange(2, self.nlayers+1):
            z = self.zs[-k]
            sp = TransferFunctionInn_derivative(z,self.TransferFunction)
            z1 = np.dot(delta, self.layer[-k+1].w.T)
            delta = z1 * sp
            self.layer[-k].db = np.sum( delta, 0 )
            self.layer[-k].dw = np.dot(self.activations[-k-1].T, delta)

        for k in range(self.nlayers):
            self.layer[k].dw += Config.C1 * self.layer[k].w - \
                                Config.C2 * np.random.rand(1) * (self.layer[k].w - self.layer[k].w_opt)  - \
                                Config.C3 * np.random.rand(1) * (self.layer[k].w - self.layer[k].w_global_opt)        
            self.layer[k].w  += Config.Weight * self.layer[k].dw       
            self.layer[k].db += Config.C1 * self.layer[k].b - \
                                Config.C2 * np.random.rand(1) * (self.layer[k].b - self.layer[k].b_opt) - \
                                Config.C3 * np.random.rand(1) * (self.layer[k].b - self.layer[k].b_global_opt)       
            self.layer[k].b  += Config.Weight * self.layer[k].db       

##########################################################################################################
    def cost_derivative(self, output_activations, y):
        """Return the vector of partial derivatives \partial C_x /
        \partial a for the output activations."""
        return (y - output_activations)

##########################################################################################################
    def rmse(self,output):
        e1 = self.cost_derivative(output, self.y)
        e  = np.sum( e1 * e1 )
        return np.sqrt( e / (e1.shape[0]*e1.shape[1]) )


           
