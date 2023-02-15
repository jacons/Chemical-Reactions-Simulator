from abc import abstractmethod


class StochasticSimulator:
    @abstractmethod
    def step(self) -> float:
        pass
