import matplotlib.pyplot as plt

from algorithms.DifferentialEq import DifferentialEq
from parsers import json_parser
from utils import solveOde

if __name__ == '__main__':

    reactions, init_state, events = json_parser(path='../sources/source2.json')
    if reactions is None or init_state is None:
        print("Runtime error in parsing phase")

    end_time = 0.001

    simulation_, times_ = solveOde(
        model_=DifferentialEq(reactions=reactions.copy(), initial_state=init_state.copy()),
        end_time=end_time,
        precision=0.00001)

    plt.grid()
    plt.title("Stochastic Chemical Reaction Simulation")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")

    plt.plot(times_, simulation_, alpha=1, linewidth=0.7)

    plt.show()
