import matplotlib.pyplot as plt

from algorithms.Tau_leaping2 import TauLeaping2
from parsers import matrix_parser
from utils import simulation, colors

if __name__ == '__main__':

    state, fields, events, mol2id = matrix_parser(path='../sources/source1.json')

    itr = 1
    end_time = 100

    hists, times = [], []
    for _ in range(itr):
        simulation_, times_ = simulation(
            model_=TauLeaping2(initial_state=state.copy(), fields=fields, mol2id=mol2id, tau=0.1),
            end_time=end_time,
            events=events.copy())

        hists.append(simulation_)
        times.append(times_)

    plt.grid()
    plt.title("Stochastic Chemical Reaction Simulation")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")

    for sim, times_ in zip(hists, times):
        for idx in range(state.shape[0]):
            plt.plot(times_, sim[:, idx], color=colors[idx], alpha=1, linewidth=0.7)
    plt.legend()
    plt.show()
