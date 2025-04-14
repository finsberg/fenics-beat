<p align="center">
  <img width="300" height="300" src="https://raw.githubusercontent.com/finsberg/fenics-beat/main/docs/_static/fenics-beat-logo.png">
</p>

[![pre-commit](https://github.com/finsberg/fenics-beat/actions/workflows/pre-commit.yml/badge.svg)](https://github.com/finsberg/fenics-beat/actions/workflows/pre-commit.yml)
[![Create and publish a Docker image](https://github.com/finsberg/fenics-beat/actions/workflows/docker-image.yml/badge.svg)](https://github.com/finsberg/fenics-beat/pkgs/container/fenics-beat)
[![Build docs](https://github.com/finsberg/fenics-beat/actions/workflows/build_docs.yml/badge.svg)](https://github.com/finsberg/fenics-beat/actions/workflows/build_docs.yml)
[![Test package](https://github.com/finsberg/fenics-beat/actions/workflows/main.yml/badge.svg)](https://github.com/finsberg/fenics-beat/actions/workflows/main.yml)
[![PyPI version](https://badge.fury.io/py/fenics-beat.svg)](https://badge.fury.io/py/fenics-beat)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13323643.svg)](https://doi.org/10.5281/zenodo.13323643)

---

# fenics-beat

Library for running cardiac electrophysiology simulations.

- Source code: https://github.com/finsberg/fenics-beat
- Documentation: https://finsberg.github.io/fenics-beat

> [!NOTE]
> This library are using the legacy version of FEnICS. New users are encourage to use [`fenicsx-beat`](https://github.com/finsberg/fenicsx-beat) which is based on FEniCSx


## Getting started

Check out the examples at https://finsberg.github.io/fenics-beat/

## Install

### Using docker (recommended)
The simplest way to use `fenics-beat` is to use the provided docker image. You can get this image by pulling it from the github registry
```
docker pull ghcr.io/finsberg/fenics-beat:latest
```
It is also possible to pull a specific version by changing the tag, e.g.
```
docker pull ghcr.io/finsberg/fenics-beat:v0.0.8
```
will use version 0.0.8.

In order to start a container you can use the [`docker run`](https://docs.docker.com/engine/reference/commandline/run/) command. For example the command
```
docker run --rm -v $(pwd):/home/shared -w /home/shared -ti ghcr.io/finsberg/fenics-beat:latest
```
will run the latest version and share your current working directory with the container.
The source code of `fenics-beat` is located at `/repo` in the docker container.

### Using pip
`fenics-beat` is also available on [pypi](https://pypi.org/project/fenics-beat/) and can be installed with
```
python3 -m pip install fenics-beat
```
However this requires FEniCS to already be installed. Currently, FEniCS can be installed by building [from source](https://bitbucket.org/fenics-project/dolfin/src/master/), using [conda](https://anaconda.org/conda-forge/fenics) or use some of the [pre-built docker images](https://github.com/orgs/scientificcomputing/packages?repo_name=packages)


## Automated tests
Upon pushing new code to the repository, a number of tests run:
* pre-commit tests.
    - Install `pre-commit`:
    ```
    python3 -m pip install pre-commit
    ```
    - Run pre-commit hooks:
    ```
    pre-commit run --all
    ```
* unit and integration tests can be found in `tests` folder
    - Install tests dependencies:
    ```
    python3 -m pip install -e .[test]
    ```
    - Run tests
    ```
    python3 -m pytest
    ```
* Examples: All examples are run as part of building the documentation

## Contributing guidelines

Detailed contributing guidelines are given [here](https://finsberg.github.io/fenics-beat/CONTRIBUTING.html).

## License
MIT

## Need help or having issues
Please submit an [issue](https://github.com/finsberg/fenics-beat/issues)
