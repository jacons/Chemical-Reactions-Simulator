import matplotlib.pyplot as plt

from algorithms.Gillespie import Gillespie
from algorithms.Tau_leaping import TauLeaping
from parsers import json_parser
from utils import simulation, avg_simulations

if __name__ == '__main__':

    reactions, init_state, events = json_parser(path='../sources/source2.json')
    if reactions is None or init_state is None:
        print("Runtime error in parsing phase")

    itr = 100
    end_time = 0.001

    hists, times = [], []
    for _ in range(itr):
        simulation_, times_ = simulation(
            model_=TauLeaping(reactions=reactions.copy(), initial_state=init_state.copy(),tau=0.01),
            end_time=end_time,
            events=events.copy())

        hists.append(simulation_)
        times.append(times_)

    for sim, times_ in zip(hists, times):
        plt.plot(times_, sim, color="gray", alpha=0.1, linewidth=1)

    x, y = avg_simulations(hists, times)
    y.sort()

    labels = [m.replace("_","") for m in init_state.keys()]
    plt.plot(y, x, marker=".", label=labels)
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")
    plt.grid()
    plt.xlim([-end_time / 100, end_time])
    plt.show()
