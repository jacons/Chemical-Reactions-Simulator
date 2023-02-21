from typing import Union

from flask import Flask, request, render_template

from Simulation import getSimulationA
from parsers import reactions_parser, state_parser, events_parser

app = Flask(__name__)

# We preload a system of chemical reaction as example
cache = {  # We save the previous computation
    "reactions": "1_S + 1_E  -> 1_SE , 0.006\n1_SE -> 1_S + 1_E , 0.005\n1_SE -> 1_P + 1_E , 0.2",
    "init_state": "E : 50 ; \nS : 200 ; \nP : 0 ;\nSE : 10 ;\nX : 0 ;",
    "events": "20 S 10 ;\n30 P -10 ;",
    "precision": 0.01,
    "end_time": 50,
    "itr": 20}


def checkinputs(form: dict) -> Union[dict, None]:
    """
    Check if the input form are filled and if the values are coherent
    Oss. events have not mandatory!
    """

    if form["reactions"] == "" or form["init_state"] == "" or form["algo"] == "" or \
            form["precision"] == "" or form["end_time"] == "" or form["itr"] == "":
        return None

    inputs = {"precision": max(0.00001, float(form["precision"])),
              "itr": max(1, int(form["itr"])),
              "end_time": max(0.0001, float(form["end_time"])),
              "reactions": form["reactions"],
              "init_state": form["init_state"],
              "events": form["events"],
              "algo": form["algo"]}

    return inputs


@app.route('/', methods=('GET', 'POST'))
def create():
    final_state = None
    if request.method == 'POST' and "run" in request.form:

        inputs = checkinputs(form=request.form)
        if inputs is not None:
            # parsing the inputs
            reactions = reactions_parser(inputs["reactions"])  # list of reactions
            init_state = state_parser(inputs["init_state"])  # dict of initial state
            events = events_parser(inputs["events"])  # list of events(tuple)

            cache["reactions"] = inputs["reactions"]
            cache["init_state"] = inputs["init_state"]
            cache["events"] = inputs["events"]
            cache["end_time"] = inputs["end_time"]
            cache["itr"] = inputs["itr"]
            cache["precision"] = inputs["precision"]

            # for graphical interface we do not take into account the matrix-based implementation
            final_state = getSimulationA(reactions, init_state, events, inputs["algo"],
                                         inputs["itr"], inputs["end_time"], inputs["precision"])

    return render_template('main.html', cache=cache, url='/static/simulation.png', final_state=final_state)
