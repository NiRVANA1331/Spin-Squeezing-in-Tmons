#module for SU(N) algebra operators
import numpy as np
from qutip import *

def ensemble_op(op,n,N):
    #returns IxIxIx...(op_n)xIxI...xI
    #I = I_(NxN)
    if n > N-1:
        raise ValueError("index needs to be between 0 and"+str(N)+"-1")

    Neye = []
    for k in  range(N):
        Neye.append(qeye(op.shape[0]))
    #(op_n)
    Neye[n] = Qobj(op)
    nid = tensor(Neye)
    return nid

def op_padded(op,N):
    #N dim I is created
    #if dim(op) = n, first nxn block on the diagonal is replaced by op
    if op.type !="oper":
        raise ValueError("entered invalid operator")
    op_N = np.identity(N, dtype=complex)
    op_dim = op.shape[0]
    op_N[0:op_dim, 0:op_dim] = op.full()
    return op_N

#create spin operators
class spin_algebra:
    #class to create all operators for a system of N spin systems
    def __init__(self, N, spin):
        #N is no of spins
        #spin is the value of spin i.e 1/2,1 ...etc
        self.N = N
        self.spin = spin
        self.dim = spin*2 
        #min number of levels in the spin
        #also used to calculate the default dimension of representation
        
    
    def generators_arr(self):
        #return all generators for a spin-n algebra
        #for now with fixed dimensions
        

        return
    
    def tensor_operator(self,op,i):
        #returns IxIxIx...(op_i)xIxI...xI
        if i > self.N-1:
            raise ValueError("index needs to be between 0 and"+str(self.N)+"-1")
        Neye = []
        for k in  range(self.N):
            Neye.append(qeye(self.spin))
        #(op_i)
        Neye[n] = basis(dim,i)*basis(dim,j).dag()
        nid = tensor(Neye)
        return nid
        
    def J_n(self,n):
        #\sum (n.J)_i
        #it should take both 3x1 and dx1 vectors
        #3x1 leads to use of just sigma x,y,z
        #dx1 leads to use of all generators
        jx = jmat(i, 'x')
        jy = jmat(i, 'y')
        jz = jmat(i, 'z')
        
        
     
    def nonlin_operators(self, degree):
        #take array of generators and combine them to get non-lin operators
        #return vector with non-lin operators only
        lin_op_arr = self.generators_arr
        nonlin_arr = []