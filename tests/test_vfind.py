from vfind import find_variants
import polars as pl


def test_vfind():
    variants = find_variants(
        "tests/test_data/toy.fq.gz",
        ("GGGCCCAGCCGGCCGGAT", "CCGGAGGCGGAGGTTCAG"),
        show_progress=False,
    )
    ground_truth = pl.read_csv("tests/test_data/ground_truth.csv")

    assert variants.equals(ground_truth)
