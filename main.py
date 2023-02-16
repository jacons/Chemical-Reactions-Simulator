from Gillespie import Gillespie
from utils import parser, single_simulation
import matplotlib.pyplot as plt

colors = {0: "b", 1: "g", 2: "r", 3: "c", 4: "m", 5: "y", 6: "k"}

if __name__ == '__main__':

    reactions, init_state, _ = parser(path='sources/source1.json')

    if reactions is None or init_state is None:
        print("Runtime error in parsing phase")

    itr = 1
    end_time = 80

    hists, times = [], []
    for _ in range(itr):
        simulation, times_ = single_simulation(
            Gillespie(reactions=reactions.copy(), initial_state=init_state.copy()), end_time)
        hists.append(simulation)
        times.append(times_)

    plt.grid()
    plt.title("Stochastic Chemical Reaction Simulation")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")

    for sim, times_ in zip(hists, times):
        for idx, m in enumerate(init_state.keys()):
            plt.plot(times_, sim[:, idx], color=colors[idx], alpha=1, linewidth=0.7)
    plt.show()
