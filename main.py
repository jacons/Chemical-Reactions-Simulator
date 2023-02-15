import numpy as np

from Gillespie import Gillespie
from utils import parser
import matplotlib.pyplot as plt

reactions, init_state = parser(path='sources/source1.json')

g = Gillespie(reactions, init_state)


def run():
    t = 0
    hist, time = [], []

    while t < 50:
        dt = g.step()
        t += dt
        time.append(t)
        hist.append(list(g.state.values()))

    hist = np.array(hist)
    for idx, m in enumerate(g.state.keys()):
        plt.plot(time, hist[:, idx], label=m)
    plt.legend()
    plt.show()


run()
# print(timeit.timeit("run()", setup="from __main__ import run", number=10)/10)
