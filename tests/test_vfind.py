from vfind import find_variants
import polars as pl
import logging
import unittest


class TestFinder(unittest.TestCase):
    """
    Test find_variants from python
    """

    def test_vfind(self):
        FORMAT = (
            "%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s"
        )
        logging.basicConfig(format=FORMAT)
        logging.getLogger().setLevel(logging.WARN)

        variants = find_variants(
            "test_data/toy.fq.gz",
            ("GGGCCCAGCCGGCCGGAT", "CCGGAGGCGGAGGTTCAG"),
            show_progress=False,
        )
        ground_truth = pl.read_csv("test_data/ground_truth.csv")

        assert variants.equals(ground_truth)


if __name__ == "__main__":
    unittest.main()
