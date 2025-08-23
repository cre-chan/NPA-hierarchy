import ncpol2sdpa as ncr
import sympy
import sys

class CHSH:
    def __init__(self):
        A =ncr.generate_operators("A",2,hermitian=True)
        B =ncr.generate_operators("B",2,hermitian=True)

        GPoly = (2 + 1/2*A[0]*B[0]+1/2*A[1]*B[0]+1/2*A[0]*B[1]-1/2*A[1]*B[1])/4
        Substitutions = {
        op**2: 1 for op in A+B
        } | {
            f*g:g*f for f in A for g in B 
        }
        self.vars = A+B 
        self.substituions = Substitutions
        self.gpoly = GPoly

    def get_relaxation(self):
        return ncr.SdpRelaxation(self.vars)

    def get_maximum_violation(self):
        return {
                "substitutions":self.substituions,
                "objective" : -self.gpoly
                }
    
    def get_check_prob(self,prob=1):
        return {
                "substitutions":self.substituions,
                "objective" : 1,
                "equalities": [self.gpoly - prob]
                }
    
class MSG:
    def __init__(self):
        F = ncr.generate_operators('F',3, hermitian=True)
        G = ncr.generate_operators('G',3,hermitian=True)
        GPoly = (9/2 + 1/2 * F[0] *G[0] + 1/2 * F[1] *G[0] \
            + 1/2 * F[2] * G[0] + 1/2 * F[0] * G[1] \
            + 1/2 * F[1] * G[1] + 1/2 * F[2] * G[1] \
            + 1/2 * F[0] * G[2] + 1/2 * F[1] * G[2] \
            + 1/2 * F[2] * G[2])/9

        Substitutions = {
        op**2: 1 for op in F+G
        } | {
            f*g:g*f for f in F for g in G 
        }
        self.vars = F+G
        self.substituions = Substitutions
        self.gpoly = GPoly

    def get_relaxation(self):
        return ncr.SdpRelaxation(self.vars)

    def get_maximum_violation(self):
        return {
                "substitutions":self.substituions,
                "objective" : -self.gpoly
                }
    
    def get_check_prob(self,prob=1):
        return {
                "substitutions":self.substituions,
                "objective" : 1,
                "equalities": [self.gpoly - prob]
                }

if __name__ == "__main__":
    mode = sys.argv[1]

    level = 3
    
    msg = MSG()

    if mode == "--solve":
        sdp = msg.get_relaxation()
        sdp.get_relaxation(level, 
                        **msg.get_maximum_violation())

        # sdp.write_to_file("MSG-1.csv")
        # sdp.write_to_file("MSG-1.dat-s")
        sdp.solve()

        print("The primial optimal value: ",sdp.primal)
    elif mode == "--check":
        sdp = msg.get_relaxation()
        sdp.get_relaxation(level, 
                        **msg.get_check_prob(1))
        sdp.solve()
        
        print("Does Alice and Bob win the game with probability 1? ", 
              sdp.status.find("optimal")!=-1)
        print(sdp.primal)
    else:
        raise ValueError("sss")