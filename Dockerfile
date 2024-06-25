FROM ghcr.io/scientificcomputing/fenics-gmsh:2024-05-30 as beat_base

ENV PYVISTA_JUPYTER_BACKEND="html"

# Requirements for pyvista
RUN apt-get update && apt-get install -y libgl1-mesa-glx libxrender1 xvfb nodejs

COPY . /repo
WORKDIR /repo

ARG TARGETPLATFORM
RUN echo "Building for $TARGETPLATFORM"
RUN if [ "$TARGETPLATFORM" = "linux/arm64" ]; then python3 -m pip install "https://github.com/finsberg/vtk-aarch64/releases/download/vtk-9.2.6-cp310/vtk-9.2.6.dev0-cp310-cp310-linux_aarch64.whl"; fi

RUN python3 -m pip install ".[demos]" jupytext


# Convert all python files to notebooks
RUN jupytext demos/*.py --to ipynb

# Jupyter-lab images for examples
FROM beat_base as beat_lab
EXPOSE 8888/tcp
ENTRYPOINT [ "jupyter", "lab", "--ip", "0.0.0.0", "--port", "8888", "--no-browser", "--allow-root"]
