import matplotlib.pyplot as plt

from algorithms.Gillespie import Gillespie
from algorithms.Tau_leaping import TauLeaping
from parsers import json_parser
from utils import simulation

colors = {0: "b", 1: "g", 2: "r", 3: "c", 4: "m", 5: "y", 6: "k"}

# qunate volte una reazion viene eseguita??

if __name__ == '__main__':

    reactions, init_state, events = json_parser(path='sources/source1.json')
    if reactions is None or init_state is None:
        print("Runtime error in parsing phase")

    itr = 1
    end_time = 100

    hists, times = [], []
    for _ in range(itr):

        simulation_, times_ = simulation(
            model_=TauLeaping(reactions=reactions.copy(), initial_state=init_state.copy(),tau=0.001),
            end_time=end_time,
            events=events.copy())

        hists.append(simulation_)
        times.append(times_)

    plt.grid()
    plt.title("Stochastic Chemical Reaction Simulation")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")

    for sim, times_ in zip(hists, times):
        for idx, m in enumerate(init_state.keys()):
            plt.plot(times_, sim[:, idx], color=colors[idx], alpha=1, linewidth=0.7)
    plt.show()
