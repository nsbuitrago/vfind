from vfind import Adapters, AlignParams, find_variants
import polars as pl
import unittest


class TestFinder(unittest.TestCase):
    """
    Test find_variants from python
    """

    def test_vfind(self):
        adapters = Adapters("GGGCCCAGCCGGCCGGAT", "CCGGAGGCGGAGGTTCAG")
        variants = find_variants("test_data/toy.fq.gz", adapters)
        ground_truth = pl.read_csv("test_data/ground_truth.csv")

        assert variants.equals(ground_truth)


if __name__ == "__main__":
    unittest.main()
