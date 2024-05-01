# Developing

vFind is written in Rust and uses [PyO3](https://pyo3.rs/v0.21.1/) and [Maturin](https://github.com/PyO3/maturin)
to generate the Python package. To get started, you will need to have the
[rust toolchain](https://www.rust-lang.org/tools/install) and [Python >= 3.10](https://www.python.org/downloads/).

Below are some general steps to get up and running. These examples
use [uv](https://github.com/astral-sh/uv). However, you could do this with
standard pip or your preferred method.

1. clone the repository to your machine

```bash
git clone git@github.com:nsbuitrago/vfind.git
cd vfind
```

2. Create a new Python virtual environment and install dev dependencies.

```bash
uv venv
source .venv/bin/activate

# this will intsall maturin, polars, and any other required packages.
uv pip install -r dev-requirements.txt
```

3. Build and install the vFind package in your virtual environment with maturin. 

```bash
# in the root project directory
maturin develop
```

4. From here, you can make changes to the Rust library and run `maturin develop`
to rebuild the package and install it in your virtual environment.

## Using the supplied flake.nix

For those that are familiar with [flakes](https://wiki.nixos.org/wiki/Flakes), there is a provided [flake.nix](/flake.nix)
that uses [devenv](https://devenv.sh) to setup a development environment with all the necessary
tools installed. To get started, make sure you have [Nix](https://nixos.org/download/) installed and have
enabled `nix-commands` and `flakes`.

To create a new shell with development tools, run `nix develop --impure`. This
may take some time on the first run since it will need to build the tools.

Once in the development shell, you can run continue with the previous steps
to make necessary changes and build the Python package.

