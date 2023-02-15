from typing import Tuple, List

import numpy as np
from numpy import ndarray

import StochasticSimulator
from Gillespie import Gillespie
from Tau_leaping import TauLeaping
from utils import parser
import matplotlib.pyplot as plt


def run_simulation(algo: StochasticSimulator, end_time: int) -> tuple[ndarray, list]:
    t = 0
    hist, times = [], []

    while t < end_time:
        dt = algo.step()
        t += dt
        times.append(t)
        hist.append(list(algo.state.values()))
    hist = np.array(hist)
    return hist, times


if __name__ == '__main__':
    reactions, init_state = parser(path='sources/source1.json')

    gillespie = Gillespie(reactions, init_state)
    simulation, time = run_simulation(gillespie, 70)

    _, ax = plt.subplots(ncols=2, figsize=(14, 5))

    for idx, m in enumerate(gillespie.state.keys()):
        ax[0].plot(time, simulation[:, idx], label=m)
    ax[0].legend()

    tau_leaping = TauLeaping(reactions, init_state, tau=0.001)
    simulation, time = run_simulation(tau_leaping, 70)

    for idx, m in enumerate(tau_leaping.state.keys()):
        ax[1].plot(time, simulation[:, idx], label=m)
    ax[1].legend()

    plt.show()
