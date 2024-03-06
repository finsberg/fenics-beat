FROM ghcr.io/scientificcomputing/fenics-gmsh:2024-02-19

COPY . /app
WORKDIR /app

RUN python3 -m pip install ".[demos]"
