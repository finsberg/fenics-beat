import abc


class CellModel:
    def __init__(
        self,
        params: dict[str, float] | None = None,
        init_conditions: dict[str, float] | None = None,
    ) -> None:
        self._parameters = type(self).default_parameters()
        if params is not None:
            self._parameters.update(params)

        self._ic = type(self).default_initial_conditions()
        if init_conditions is not None:
            self._ic.update(init_conditions)

    @property
    def ic(self) -> dict[str, float]:
        """Initial conditions"""
        return self._ic

    @property
    def parameters(self) -> dict[str, float]:
        """Model parameters"""
        return self._parameters

    @staticmethod
    @abc.abstractmethod
    def default_parameters() -> dict[str, float]:
        ...

    @staticmethod
    @abc.abstractmethod
    def default_initial_conditions() -> dict[str, float]:
        ...

    @abc.abstractmethod
    def F(self, v, s, time=None):
        "Return right-hand side for state variable evolution."
        ...

    @abc.abstractmethod
    def I(self, v, s, time=None):
        "Return the ionic current."
        ...

    @abc.abstractmethod
    def num_states(self):
        """Return number of state variables (in addition to the
        membrane potential)."""
        ...
