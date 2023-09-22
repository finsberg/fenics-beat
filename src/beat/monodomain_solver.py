import logging
from dataclasses import dataclass
from .monodomain_model import MonodomainModel
from .odesolver import DolfinODESolver

logger = logging.getLogger(__name__)
EPS = 1e-12


@dataclass
class MonodomainSplittingSolver:
    pde: MonodomainModel
    ode: DolfinODESolver

    def solve(self, interval, dt):
        (T0, T) = interval
        if dt is None:
            dt = T - T0
        t0 = T0
        t1 = T0 + dt

        while t1 < T + EPS:
            logger.info(f"Solving on t = ({t0:.2f}, {t0:.2f})")
            self.step((t0, t1))

            t0 = t1
            t1 = t0 + dt

    def step(self, interval):
        # Extract some parameters for readability
        theta = self.pde.parameters["theta"]

        # Extract time domain
        (t0, t1) = interval
        dt = t1 - t0
        t = t0 + theta * dt

        logger.info("Tentative ODE step")

        # Solve ODE
        self.ode.step(t0, theta * dt)
        # Move voltage to FEniCS
        self.ode.to_dolfin()

        logger.info("PDE step")
        # Solve PDE
        self.pde.step((t0, t1))

        # Copy voltage from PDE to ODE
        self.ode.from_dolfin()
        # If first order splitting, we are done.
        if theta == 1.0:
            return

        # Otherwise, we do another ode_step:
        logger.info("Corrective ODE step")

        # To the correction step
        self.ode.step(t, (1.0 - theta) * dt)
        # And copy the solution back to FEniCS
        self.ode.to_dolfin()
