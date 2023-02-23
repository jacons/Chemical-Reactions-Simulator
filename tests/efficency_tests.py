import time

from Simulation import getSimulationA, getSimulationB
from parsers import dicts_parser, matrix_parser


def executions(path='../sources/source1.json', itrs=range(10, 100, 5), end_time=100):
    reactions, init_state, events = dicts_parser(path=path)
    state, fields, events, mol2id = matrix_parser(path=path)

    for itr in itrs:
        times = []
        print(itr, "\t", end="")

        start_time = time.time()
        r = getSimulationA(reactions, init_state, events, "Gillespie", itr, end_time, 0, figure=False)
        t = (time.time() - start_time)
        times.append(t)

        start_time = time.time()
        r = getSimulationB(state, fields, events, mol2id, "Gillespie2", itr, end_time, 0)
        t = (time.time() - start_time)
        times.append(t)

        start_time = time.time()
        r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.1, figure=False)
        t = (time.time() - start_time)
        times.append(t)

        start_time = time.time()
        r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.01, figure=False)
        t = (time.time() - start_time)
        times.append(t)

        start_time = time.time()
        r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.1)
        t = (time.time() - start_time)
        times.append(t)

        for t in times:
            print(round(t, 4), "\t", end="")
        print()


def run_test(path='../sources/source1.json'):
    itr = 10
    end_time = 40
    reactions, init_state, events = dicts_parser(path=path)
    state, fields, events, mol2id = matrix_parser(path=path)

    # --------------------
    print("Gillespie method")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "Gillespie", itr, end_time, 0, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.5")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.5, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.1")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.1, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.01")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.01, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.001")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.001, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("ODE with 0.1")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision=0.1, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("ODE with 0.01")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision=0.01, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("ODE with 0.001")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision=0.001, figure=False)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("Gillespie matrix based")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "Gillespie2", itr, end_time, 0)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.5")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.5)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.1")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.1)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.01")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.01)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.001")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.001)
    t = (time.time() - start_time)
    print("Ends in --- ", "{:.2f}".format(t), "seconds --- ", r, "\n")
    # --------------------


run_test()
