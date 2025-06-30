from vfind import find_variants, align
import polars as pl
from polars.testing import assert_frame_equal
import os.path as path

TEST_DATA_PATH = path.join("tests", "test_data")

MAGICAL_DF = pl.DataFrame([
    pl.Series("sequence", ["MAGICAL"], dtype=pl.String),
    pl.Series("count", [4], dtype=pl.UInt64)
])

MAGICAL_ADAPTERS_DF = pl.DataFrame([
    pl.Series("sequence", [
        "GPSRRDMAGICALPEAEVQ",
        "GPSRRDMAGICALPEAEVE",
        "GPSRPDMAGICALPEAEVQ",
        "GPSRPDMAGICALPEAEVE",
    ], dtype=pl.String),
    pl.Series("count", [1,1,1,1], dtype=pl.UInt64)

])

LIIGAND_DF = pl.DataFrame([
    pl.Series("sequence", ["LIIGAND"], dtype=pl.String),
    pl.Series("count", [4], dtype=pl.UInt64)
])

def _test_recovery(adapters: tuple[str, str], fq_path: str, expected_df: pl.DataFrame) -> None:
    """
    Helper function for testing variant recovery

    Parameters
    ----------

    adapters (tuple(str, str)): Prefix and suffix adapters to search.
    fq_path (str): Path to input fastq file.
    csv_path (str): Path to ground truth csv file.
    """
    variants = find_variants(fq_path, adapters, show_progress=False)
    assert_frame_equal(variants, expected_df, check_row_order=False)
    #assert variants.equals(expected_df)


def test_long_adapter_recovery():
    _test_recovery(
        ("GGGCCCAGCCGGCCGGAT", "CCGGAGGCGGAGGTTCAG"),
        path.join(TEST_DATA_PATH, "toy_18bp_barcode.fq.gz"),
        #path.join(TEST_DATA_PATH, "toy_ground_truth.csv")
        MAGICAL_DF
    )


def test_short_adapter_recovery():
    _test_recovery(
        ("GATCATG", "GAACTGC"),
        path.join(TEST_DATA_PATH, "toy_7bp_barcode.fq.gz"),
        MAGICAL_DF
    )


def test_magical_demultiplex():
    _test_recovery(
        ("GATCATG", "GAACTGC"),
        path.join(TEST_DATA_PATH, "demultiplex.fq.gz"),
        MAGICAL_DF
    )

    _test_recovery(
        ("GATCATG", "ACCAGGT"),
        path.join(TEST_DATA_PATH, "demultiplex.fq.gz"),
        LIIGAND_DF
    )


def test_keep_adapters():
    variants = find_variants(
        path.join(TEST_DATA_PATH, "toy_18bp_barcode.fq.gz"),
        ("GGGCCCAGCCGGCCGGAT", "CCGGAGGCGGAGGTTCAG"),
        skip_trimming=True,
        show_progress=False,
    )

    assert_frame_equal(variants, MAGICAL_ADAPTERS_DF, check_row_order=False)
