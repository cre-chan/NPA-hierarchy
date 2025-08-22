import sdpap
import sys


usage ="""
python run-sdpa.py path_to_input_file path_to_output_file
"""

class SolverArgs:
    def __init__(self):
        prob_filename, result_filename = sys.argv[1:3]
        self.infile = prob_filename
        self.outfile = result_filename
        print(prob_filename,result_filename)

    @property
    def sdpopt(self):
        return {
            "print":"no",
            "sdpaResult":self.outfile
        }


if __name__ == "__main__":
    args = SolverArgs()
    filename = args.infile
    A, b, c, K, J = sdpap.importsdpa(filename)
    x, y, sdpapinfo, timeinfo, sdpainfo = sdpap.solve(A,b,c,K,J,args.sdpopt)
