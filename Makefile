dev:
	uv run maturin develop --uv
tests:
	echo "Running cargo unit tests..." \
	&& cargo test \
	&& echo "Running python integration tests..." \
	&& uv run maturin develop --uv \
	&& uv run python3 -m unittest tests/test_vfind.py
