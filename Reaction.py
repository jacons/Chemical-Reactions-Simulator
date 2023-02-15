class Reaction:
    def __init__(self, name: str, reactants: dict, products: dict, kinetic: float):
        self.name = name

        self.reactants = reactants
        self.products = products

        self.kinetic = kinetic

    def show_reaction(self) -> str:
        reaction = ""
        for r, l in self.reactants.items():
            reaction += str(l) + str(r) + " "
        reaction += " --> "
        for r, l in self.products.items():
            reaction += str(l) + str(r) + " "

        return reaction + " (" + str(self.kinetic) + ")"