# Import packages.
import cvxpy as cp
import numpy as np
import scipy.io as sio
import scipy.sparse as ssparse
import h5py

# The file containing the description of the SDP programs
inputFilePath = "./CHSH/CDMSG-NPA-2.mat" 

def toSparseMat(g,dim):
    data = g["data"][:]
    ir = g["ir"][:]
    jc = g["jc"][:]

    return ssparse.csc_matrix((data, ir, jc), shape=dim)

def readMat(path):
    try:
        file = sio.loadmat(path)
        PG = file["PG"]
        UL = file["UL"]
        Cons = file["CK"].tocsc()
    except NotImplementedError:
        with h5py.File(path,"r") as file:
            width = file["width"][:].item()
            numCons = file["numCons"][:].item()
            PG = toSparseMat(file["PG"],(width,width))
            UL = toSparseMat(file["UL"],(width,width))
            Cons = toSparseMat(file["CK"],(numCons*width,width))
    
    return (PG,UL,Cons)

def readSDP(path):
    print("Reading file ",path)
    PG,UL,Cons = readMat(path)
    print("Succeeded in reading file")
    h,w = Cons.shape
    print("The dimension of CK is ",Cons.shape)
    CK = []
    for i in range(h//w):
        start = i*w
        end = (i+1)*w 
        CK.append(Cons[start:end,:])
    
    return (PG,UL,CK)

# CHSH = sio.loadmat(inputFilePath)

# PG is the coefficent matrix
# PG = CHSH["PG"]
PG,UL,CKs = readSDP(inputFilePath)

print("The dimension of PG is",PG.shape)
print("The dimension of CK is",CKs[0].shape)
NPAvars = cp.Variable(
    PG.shape,
    symmetric = True
)

objFunc = cp.trace(
    PG @ NPAvars
)

# consList = CHSH["ConsKeys"]
Cons = [NPAvars >> 0]
Cons += [
    cp.trace(mat @ NPAvars) == 0 for mat in CKs
]
Cons +=[
    cp.trace(UL @ NPAvars) == 1
]

primal = cp.Problem(
    cp.Maximize(objFunc),
    Cons
)

print("Solving problem")
primal.solve()
print("Succeeded")

print("The optimal value is ", primal.value)
# print("A solution X is")
# print(NPAvars.value)