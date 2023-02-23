import time

import matplotlib.pyplot as plt

from algorithms.Gillespie import Gillespie
from parsers import dicts_parser
from utils import simulation, avg_simulations

if __name__ == '__main__':

    reactions, init_state, events = dicts_parser(path='../sources/source1.json')
    if reactions is None or init_state is None:
        print("Runtime error in parsing phase")

    itr = 100
    end_time = 50

    hists, times = [], []
    start_time = time.time()

    for _ in range(itr):
        simulation_, times_ = simulation(
            model_=Gillespie(reactions=reactions.copy(), initial_state=init_state.copy()),
            end_time=end_time,
            events=events.copy())

        hists.append(simulation_)
        times.append(times_)

    t = (time.time() - start_time)
    print("Ends in --- ", "{:.4f}".format(t), "seconds --- ", "\n")

    for sim, times_ in zip(hists, times):
        plt.plot(times_, sim, color="gray", alpha=0.1, linewidth=1)

    x, y = avg_simulations(hists, times)
    y.sort()

    labels = [m.replace("_", "") for m in init_state.keys()]
    plt.plot(y, x, label=labels)
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")
    plt.grid()
    plt.xlim([-end_time / 100, end_time])
    plt.savefig("fig.png")
    plt.show()
