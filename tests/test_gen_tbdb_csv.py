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
    ("embA_2045_del_12_accgcatcttggc_a", "embA", "c.2046_2057del"),  # deletion of multiple bases on the positive strand
    ("gid_103_del_1_gc_g", "gid", "c.102del"),  # deletion of a single base in the negative strand
    ("pncA_390_del_4_cacat_c", "pncA", "c.386_389del"),  # deletion of multiple bases on the negative strand
    ("rpoB_1296_ins_3_a_attc", "rpoB", "c.1296_1297insTTC"),  # insertion of multiple bases in the positive strand
    ("pncA_392_ins_1_a_ac", "pncA", "c.391_392insC"),  # insertion of a single base on the negative strand with offset
    ("pncA_533_ins_1_gcggtgcgcatctcctccagcgcggcgacggtgg_gcggtgcgcatctcctcccagcgcggcgacggtgg", "pncA", "c.515_516insC"), # insertion of a single base on the negative strand with offset
    ("pncA_-4_del_1_tc_t", "pncA", "c.-5del"),  # promoter deletion
    ("whiB6_-94_del_3_gagt_g", "whiB6", "c.-95_-97del"),  # promoter deletion of multiple bases
    ("whiB6_-66_del_1_agctccgagctctagt_agctccgagcctagt", "whiB6", "c.-76del"),  # promoter deletion with an offset
    ("gyrB_-60_ins_1_t_tg", "gyrB", "c.-60_-61insG"),  # promoter insertion single base with no offset
    ("gyrB_-165_ins_2_cg_cgtc", "gyrB", "c.-166_-167insTC"),  # promoter insertion of two bases with an offset of 1
    ("rrs_-26_ins_10_tcccttt_tcccttttccaaaggga", "rrs", "c.-32_-33insTCCAAAGGGA"),  # promoter insertion with offset of 7 and multiple bases
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

