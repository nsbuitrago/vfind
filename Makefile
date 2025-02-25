dev:
	@uv run maturin dev --uv
test:
	@echo "--- Running rust unit tests ---" && \
	cargo test --verbose && \
	echo "--- Running python integration tests ---" && \
	uv run pytest
