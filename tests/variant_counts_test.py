import filecmp
import os
import pytest
import subprocess

from triphecta import variant_counts


def test_save_and_load_tsv_files():
    tmp_file = "tmp.variant_counts.tsv"
    subprocess.check_output(f"rm -f {tmp_file}", shell=True)
    variants = [
        variant_counts.VariantCounts(het=1, hom=2, null=3, het_to_hom=4),
        variant_counts.VariantCounts(het=42, hom=43, null=44, het_to_hom=45),
    ]
    variant_counts.save_variant_count_list_to_tsv(variants, tmp_file)
    got_variants = variant_counts.load_variant_count_list_from_tsv(tmp_file)
    assert got_variants == variants
    os.unlink(tmp_file)
