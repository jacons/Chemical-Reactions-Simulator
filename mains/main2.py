import matplotlib.pyplot as plt

from algorithms.Tau_leaping2 import TauLeaping2
from parsers import matrix_parser
from utils import simulation, avg_simulations

if __name__ == '__main__':

    init_state, fields, events, mol2id = matrix_parser(path='../sources/source1.json')

    itr = 100
    end_time = 50

    hists, times = [], []
    for _ in range(itr):
        simulation_, times_ = simulation(
            model_=TauLeaping2(initial_state=init_state.copy(), fields=fields, mol2id=mol2id, tau=0.5),
            end_time=end_time,
            events=events.copy())

        hists.append(simulation_)
        times.append(times_)

    for sim, times_ in zip(hists, times):
        plt.plot(times_, sim, color="gray", alpha=0.1, linewidth=1)

    x, y = avg_simulations(hists, times)
    y.sort()

    labels = [m.replace("_", "") for m in mol2id.keys()]
    plt.plot(y, x, label=labels)
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("N° Molecules")
    plt.grid()
    plt.xlim([-end_time / 100, end_time])
    plt.show()
