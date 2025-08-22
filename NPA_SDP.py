# Import packages.
import cvxpy as cp
import numpy as np
import scipy.io as sio
import scipy.sparse as ssparse
import h5py
import sys
import sdpap

# The file containing the description of the SDP programs



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
        try:
            Cons = file["CK"].tocsc()
        except AttributeError:
            Cons = file["CK"]
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
    CK = []
    for i in range(h//w):
        start = i*w
        end = (i+1)*w 
        CK.append(Cons[start:end,:])
    
    return (PG,UL,CK)

class SolverArgs:
    def __init__(self):
        self.mode = "optimize"

        args = sys.argv[1:]
        it = iter(args)
        self.__get_infile(it)
        while True:
            try:
                arg = next(it)

                if not self.is_args(arg):
                    raise ValueError(f"{arg} is not a valid argument")
                attr = self.get_argname(arg)

                try:
                    func = getattr(self,f"__get_{attr}")
                except AttributeError:
                    raise ValueError(f"{arg} is not a valid argument") from None
                
                func(it)
                
            except StopIteration:
                break

    def is_args(self,str):
        return str.startswith("--")
    
    def get_argname(self,str):
        return str.lstrip("-")

    def __get_infile(self,it):
        self.infile = next(it)

    def __get_mode(self,it):
        v = next(it)
        if v=="optimize":
            # this mode finds the optimal winning probability by NPA hierarchy
            self.mode = "opt"
        elif v=="check":
            self.mode = "check"
        else:
            raise ValueError("--mode must be followed by one of optimize or check")

    def get_options(self):
        return {
            "infile":self.infile,
            "mode":self.mode
        }

if __name__ == "__main__":
    solver_args = SolverArgs()
    opts = solver_args.get_options()
    inputFilePath = opts["infile"]
    PG,UL,CKs = readSDP(inputFilePath)

    NPAvars = cp.Variable(
        PG.shape,
        symmetric = True
    )

    objFunc = cp.trace(
        PG @ NPAvars
    )

    Cons = [NPAvars >> 0]
    Cons += [
        cp.trace(mat @ NPAvars) == 0 for mat in CKs
    ]
    Cons +=[
        cp.trace(UL @ NPAvars) == 1
    ]

    print("The dimension of NPAvars: ",PG.shape)
    print("The number of constraints: ",len(CKs)+1)

    primal = cp.Problem(
        cp.Maximize(objFunc),
        Cons
    )

    print("Solving problem")
    primal.solve(solver=cp.SCS, eps=1e-8)
    print("Succeeded")

    print("The optimal value is ", primal.value)
    