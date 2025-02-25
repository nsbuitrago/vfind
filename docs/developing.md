# Developing

vFind is written in Rust and uses [PyO3](https://pyo3.rs/v0.21.1/) and
[Maturin](https://github.com/PyO3/maturin) to generate the Python package. To
get started, you will need to have the [rust
toolchain](https://www.rust-lang.org/tools/install) >= 1.83 and [Python >=
3.10](https://www.python.org/downloads/). Note that you might be able to use an
older toolchain < 1.83 but we have not tested this specifically.

Below are some general steps to get up and running. These examples
use [uv](https://github.com/astral-sh/uv). However, you could do this with
standard pip or your preferred method.

1. clone the repository to your machine

```bash
git clone git@github.com:nsbuitrago/vfind.git
cd vfind
```

2. Sync dependencies with uv and build/install vFind package

```bash
uv sync
uv run maturin develop --uv

# or with the provided Makefile
make dev
```

3. From here, you can make changes to the Rust library and rebuild/install
the package for testing on the python side.

## Using the supplied flake.nix

For those that are familiar with [flakes](https://wiki.nixos.org/wiki/Flakes),
there is a provided [flake.nix](/flake.nix) to help setup a development
environment with all the necessary tools installed. To get started, make sure
you have Nix installed and have enabled `nix-commands` and `flakes`. We
recommend using the Nix installer from [Determinate
systems](https://determinate.systems/nix-installer/) to handle this setup

To create a new shell, run `nix develop` in the root of the project. This may
take some time on the first run since it will need to build the listed packages.

Once in the development shell, you can run continue with the previous steps to
make necessary changes and build the Python package.
