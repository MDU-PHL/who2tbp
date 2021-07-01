"""
Test conversion of var from WHO to style to HGVS
"""

import pytest
from who2tbp.lib.gen_tbdb_csv import who2hgvs
from who2tbp.lib.gen_tbdb_csv import insertion_calc


@pytest.mark.parametrize("var,gene,hgvs", [
    ("rpoB_S450L", "rpoB", "p.Ser450Leu"),
    ("gid_Q125!", "gid", "p.Gln125*"),
    ("rrs_g1484t", "rrs", "r.1484g>t"),
    ("eis_g-37t", "eis", "c.-37G>T"),
    ("gid_103_del_1_gc_g", "gid", "c.104del"),
    ("pncA_390_del_4_cacat_c", "pncA", "c.391_394del"),
    ("pncA_392_ins_1_a_ac", "pncA", "c.392_393insC")
])
def test_who2hgvs(var, gene, hgvs):
    obs_gene, obs_hgvs = who2hgvs(var)
    assert obs_gene == gene
    assert obs_hgvs == hgvs


@pytest.mark.parametrize("ref,alt,ins,start_pos", [
    ("gcggtgcgcatctcctccagcgcggcgacggtgg", "gcggtgcgcatctcctcccagcgcggcgacggtgg", "C", 18),
    ("agcaccctggtggccaa", "agcaccctggtgggccaa", "G", 13),
    ("a", "acc", "CC", 1),
    ("a", "ac", "C", 1)
])
def test_insertion_calc(ref, alt, ins, start_pos):
    obs_start_pos, obs_ins = insertion_calc(ref, alt)
    assert obs_ins == ins
    assert obs_start_pos == start_pos

