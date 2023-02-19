from flask import Flask, request, render_template

from Simulation import getSimulationA
from parsers import reactions_parser, state_parser, events_parser

app = Flask(__name__)

cache = {
    "reactions": "2_A + 1_B -> 2_C , 5 \n2_C -> 2_A + 1_B , 0.5",
    "init_state": "A : 100 ; \nB : 100 ; \nC : 0 ;",
    "events": ""}


@app.route('/', methods=('GET', 'POST'))
def create():

    if request.method == 'POST':

        if "run" in request.form:

            str_reactions = request.form["reactions"]
            str_state = request.form["init_state"]
            str_events = request.form["events"]

            cache["reactions"] = str_reactions
            cache["init_state"] = str_state
            cache["events"] = str_events

            reactions = reactions_parser(str_reactions)
            init_state = state_parser(str_state)
            events = events_parser(str_events)

            algo = request.form["algo"]

            precision = float(request.form["precision"])
            end_time = float(request.form["end_time"])
            itr = int(request.form["itr"])

            getSimulationA(reactions, init_state, events, algo, itr, end_time, precision)

    return render_template('main.html', cache=cache, url='/static/simulation.png')
