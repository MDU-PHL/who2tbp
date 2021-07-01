"""
Create a TBProfiler CSV file ready for use with parse_db.py

For insertions and deletions we need to account for the strand when calculating
the position.
"""

import re
import difflib
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.Data.IUPACData import protein_letters
from typing import List


# defining regex patterns for each variant case
# note the WHO is using ! to symbolise stop codon rather than *
AA_CHANGE = re.compile(f'([{protein_letters}])([0-9]{{1,4}})([{protein_letters}!])')
ncRNA_CHANGE = re.compile(f'([acgt])([0-9]{{1,4}})([acgt])')
PROMOTER_CHANGE = re.compile(f"([acgt])(-[0-9]{{1,3}})([acgt])")
NUC_DELETION = re.compile(f"([0-9]{{1,4}})_(del)_([0-9]{{1,2}})_([actg]+)_([actg]+)")
NUC_INSERTION = re.compile(f"([0-9]{{1,4}})_(ins)_([0-9]{{1,3}})_([acgt]+)_([acgt]+)")


def insertion_calc(ref_seq: str, alt_seq:str) -> str:
    """
    Given a reference sequence and an alternate sequence, identify
    the inserted bases
    :param ref_seq: the sequence in the reference genome
    :param alt_seq: the sequence in the mutated genome
    :return: inserted sequences
    """

    diff = difflib.ndiff(ref_seq.upper(), alt_seq.upper())
    insert_data = [(pos, base.strip("+").strip()) for pos, base in enumerate(diff) if base.startswith('+')]
    insert_seq = ''.join([rec[1] for rec in insert_data])
    start_pos = insert_data[0][0]
    return start_pos, insert_seq


def who2hgvs(var: str) -> tuple:
    """
    Given a variant in the format from the WHO, return the equivalent variant
    in the HGVS nomenclature

    Some detailed exmaples from the TBDB github README:
    Amino acid substitutions. Example: S450L in rpoB would be p.Ser450Leu
    Deletions in genes. Example: Deletion of nucleotide 758 in tlyA would be c.758del
    Insertion in genes. Example: Insertion of GT between nucleotide 1850 and 1851 in katG would be c.1850_1851insGT
    SNPs in non-coding RNAs. Example: A to G at position 1401 in rrs would be r.1401a>g
    SNPs in gene promoters. Example: A to G 7 bases 5' of the start codon in pncA c.-7A>G

    :param var: in the WHO format
    :return: gene ID, and the var in the HGVS format
    """
    gene, var = var.split("_", 1)

    match_aa = AA_CHANGE.match(var)
    if match_aa:
        ref_aa = protein_letters_1to3[match_aa.group(1)]
        pos = match_aa.group(2)
        alt_aa = protein_letters_1to3[match_aa.group(3)] if match_aa.group(3) != '!' else "*"
        return gene, f"p.{ref_aa}{pos}{alt_aa}"

    match_ncrna = ncRNA_CHANGE.match(var)
    if match_ncrna:
        ref_nuc = match_ncrna.group(1)
        pos = match_ncrna.group(2)
        alt_nuc = match_ncrna.group(3)
        return gene, f"r.{pos}{ref_nuc}>{alt_nuc}"

    match_prom = PROMOTER_CHANGE.match(var)
    if match_prom:
        ref_nuc = match_prom.group(1)
        pos = match_prom.group(2)
        alt_nuc = match_prom.group(3)
        return gene, f"c.{pos}{ref_nuc.upper()}>{alt_nuc.upper()}"

    # deletions are given upstream from the indicated position
    # so, need to massage the output
    match_del = NUC_DELETION.match(var)
    if match_del:
        start_pos = int(match_del.group(1)) + 1
        end_pos = int(match_del.group(3)) + start_pos - 1
        if end_pos - start_pos == 0:
            return gene, f"c.{start_pos}del"
        else:
            return gene, f"c.{start_pos}_{end_pos}del"

    # insertions are a bit weird because i need to figure out what is the
    # base diff as the WHO report why more than required

    match_ins = NUC_INSERTION.match(var)
    if match_ins:
        indicator_pos = int(match_ins.group(1))
        ref_seq = match_ins.group(4)
        alt_seq = match_ins.group(5)
        ins_pos, ins_seq = insertion_calc(ref_seq, alt_seq)
        start_pos = indicator_pos + ins_pos
        end_pos = start_pos + len(ins_seq)
        return gene, f"c.{start_pos}_{end_pos}ins{ins_seq}"


def who2tbdw(data: List[dict]) -> dict:
    """
    Given the output from WHO extraction, translate the variant code to 
    something the HGVS nomenclature

    :param data: 
    :return: 
    """
