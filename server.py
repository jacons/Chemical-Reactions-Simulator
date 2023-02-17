from flask import Flask, request, render_template

from parsers import reactions_parser

app = Flask(__name__)


@app.route('/', methods=('GET', 'POST'))
def create():
    if request.method == 'POST':

        if "run" in request.form:
            reactions = request.form["reactions"]
            reactions = reactions_parser(reactions)

    return render_template('main.html')
