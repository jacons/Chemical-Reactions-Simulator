import time

from Simulation import getSimulationA, getSimulationB
from parsers import json_parser, matrix_parser


def run_test(path='../sources/source1.json'):
    itr = 10
    end_time = 40

    reactions, init_state, events = json_parser(path=path)
    state, fields, events, mol2id = matrix_parser(path=path)

    # --------------------
    print("Gillespie method")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "Gillespie", itr, end_time, 0, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.5")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.5, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.1")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.1, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.01")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.01, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping with 0.001")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, 0.001, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("ODE with 0.1")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision=0.1, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("ODE with 0.01")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision=0.01, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("ODE with 0.001")
    start_time = time.time()
    r = getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision=0.001, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("Gillespie matrix based")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "Gillespie2", itr, end_time, 0, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")

    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.5")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.5, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.1")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.1, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.01")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.01, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------

    # --------------------
    print("TauLeaping matrix based precision 0.001")
    start_time = time.time()
    r = getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, 0.001, figure=False)
    print("Ends in --- ", "{:.2f}".format((time.time() - start_time)), "seconds --- ", r, "\n")
    # --------------------


run_test()
