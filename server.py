from flask import Flask, request, render_template

from Simulation import getSimulationA
from parsers import reactions_parser, state_parser, events_parser

app = Flask(__name__)

cache = {  # We save the previous computation
    "reactions": "1_S + 1_E  -> 1_SE , 0.006\n1_SE -> 1_S + 1_E , 0.005\n1_SE -> 1_P + 1_E , 0.2",
    "init_state": "E : 50 ; \nS : 200 ; \nP : 0 ;\nSE : 10 ;\nX : 0 ;",
    "events": "",
    "precision": 0.1,
    "end_time": 80,
    "itr": 1}


def checkinputs(form):
    if form["reactions"] == "" or form["init_state"] == "" or form["algo"] == "" or \
            form["precision"] == "" or form["end_time"] == "" or form["itr"] == "":
        return False
    else:
        return True


@app.route('/', methods=('GET', 'POST'))
def create():
    final_state = None
    if request.method == 'POST':
        if "run" in request.form or checkinputs(request.form):
            str_reactions = request.form["reactions"]
            str_state = request.form["init_state"]
            str_events = request.form["events"]
            algo = request.form["algo"]

            precision = float(request.form["precision"])
            end_time = float(request.form["end_time"])
            itr = int(request.form["itr"])

            error_parser = False
            try:
                reactions = reactions_parser(str_reactions)
                init_state = state_parser(str_state)
                events = events_parser(str_events)
            except:
                error_parser = True

            if not error_parser:
                cache["reactions"] = str_reactions
                cache["init_state"] = str_state
                cache["events"] = str_events
                cache["end_time"] = end_time
                cache["itr"] = itr
                cache["precision"] = precision

                final_state = getSimulationA(reactions, init_state, events, algo, itr, end_time, precision)

    return render_template('main.html', cache=cache, url='/static/simulation.png', final_state=final_state)
