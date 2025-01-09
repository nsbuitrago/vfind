from vfind import find_variants
import polars as pl
import unittest


class TestFinder(unittest.TestCase):
    """
    Test find_variants from python
    """

    def test_vfind(self):
        variants = find_variants(
            "tests/data/toy.fq.gz",
            ("GGGCCCAGCCGGCCGGAT", "CCGGAGGCGGAGGTTCAG"),
            show_progress=False,
        )
        ground_truth = pl.read_csv("tests/data/ground-truth.csv")
        print(variants.head())

        assert variants.equals(ground_truth)


if __name__ == "__main__":
    unittest.main()
