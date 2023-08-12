#module for SU(N) algebra operators
import numpy as np
from qutip import *

def tensor_op(op,n,N):
    #returns IxIxIx...(op_n)xIxI...xI
    #I has dimensons of op
    if n > N-1:
        raise ValueError("index needs to be between 0 and"+str(N)+"-1")

    Neye = []
    for k in  range(N):
        Neye.append(qeye(op.shape[0]))
    #(op_n)
    Neye[n] = Qobj(op)
    nid = tensor(Neye)
    return nid

def ensemble_op(op, N):
    #sum of tensor operators
    op_ensemble = 0
    for i in range(N):
        op_ensemble += tensor_op(op,i,N)
    return op_ensemble


def op_padded(op,N):
    #N dim I is created
    #if dim(op) = n < N, first nxn block on the diagonal is replaced by op
    if op.type !="oper":
        raise ValueError("entered invalid operator")
    op_dim = op.shape[0]
    if op_dim > N:
        raise ValueError("operator dimension larger than dimension provided")
    op_N = np.identity(N, dtype=complex)
    op_N[0:op_dim, 0:op_dim] = op.full()
    return op_N

def spin1_op(op_type = "gellmann"):
    if op_type == "gellmann":
        return
    return
#dumb way for spin 1
def gellmann_matrices():
    # Gell-Mann matrices in 3x3 dimension
    lambda1 = Qobj(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex))
    lambda2 = Qobj(np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex))
    lambda3 = Qobj(np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex))
    lambda4 = Qobj(np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex))
    lambda5 = Qobj(np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex))
    lambda6 = Qobj(np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex))
    lambda7 = Qobj(np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex))
    lambda8 = Qobj(np.array([[1/np.sqrt(3), 0, 0], [0, 1/np.sqrt(3), 0], [0, 0, -2/np.sqrt(3)]], dtype=complex))

    gellmann_matrices_list = [lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8]

    return gellmann_matrices_list

def spin1_tensorop():
    #tensor operators of spin 1 as explained in Begzjav Et.al, 2021 
    sx = Qobj(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex))
    sy = Qobj(np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex))
    sz = Qobj(np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex))
    q_x2_y2 = sx**2 - sy**2
    q_3z2_r2 = 3*sz**2-2* np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=complex)
    qxy = sx*sy + sy*sx
    qyz = sy*sz + sz*sx 
    qzx = sz*sx + sx*sz
    return [sx,sy,sz,q_x2_y2, q_3z2_r2,qxy,qyz,qzx]

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