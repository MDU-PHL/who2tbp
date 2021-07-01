"""
Test conversion of var from WHO to style to HGVS
"""

import pytest
from who2tbp.lib.gen_tbdb_csv import who2hgvs
from who2tbp.lib.gen_tbdb_csv import insertion_calc
from who2tbp.lib.gen_tbdb_csv import get_gene_strands


gene_strands = get_gene_strands()


@pytest.mark.parametrize("var,gene,hgvs", [
    ("rpoB_S450L", "rpoB", "p.Ser450Leu"),
    ("gid_Q125!", "gid", "p.Gln125*"),
    ("rrs_g1484t", "rrs", "r.1484g>t"),
    ("eis_g-37t", "eis", "c.-37G>T"),
    ("embA_2045_del_12_accgcatcttggc_a", "embA", "c.2046_2057del"),  # positive strand
    ("gid_103_del_1_gc_g", "gid", "c.102del"),  # negative strand
    ("pncA_390_del_4_cacat_c", "pncA", "c.386_389del"),  # negative strand
    ("rpoB_1296_ins_3_a_attc", "rpoB", "c.1296_1297insTTC"),  # positive strand
    ("pncA_392_ins_1_a_ac", "pncA", "c.391_392insC"),  # negative strand with offset
    ("pncA_533_ins_1_gcggtgcgcatctcctccagcgcggcgacggtgg_gcggtgcgcatctcctcccagcgcggcgacggtgg", "pncA", "c.515_516insC")  # negative strand with offset
])
def test_who2hgvs(var, gene, hgvs):
    obs_gene, obs_hgvs = who2hgvs(var, gene_strands)
    assert obs_gene == gene
    assert obs_hgvs == hgvs


@pytest.mark.parametrize("ref,alt,ins,offset", [
    ("gcggtgcgcatctcctccagcgcggcgacggtgg", "gcggtgcgcatctcctcccagcgcggcgacggtgg", "C", 17),
    ("agcaccctggtggccaa", "agcaccctggtgggccaa", "G", 12),
    ("a", "acc", "CC", 0),
    ("a", "ac", "C", 0)
])
def test_insert_calc(ref, alt, ins, offset):
    obs_offset, obs_ins = insertion_calc(ref, alt)
    assert obs_ins == ins
    assert obs_offset == offset


def test_get_gene_strands():
    genes = get_gene_strands()
    assert len(genes) == 1951

