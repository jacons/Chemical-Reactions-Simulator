from parsers import json_parser


class DifferentialEq:
    def __init__(self, reactions: list, initial_state: dict):
        self.reactions: list = reactions
        self.state: dict = initial_state


reactions, init_state, events = json_parser(path='sources/source1.json')
ode = DifferentialEq(reactions, init_state)
