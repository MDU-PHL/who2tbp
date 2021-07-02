"""
Create a TBProfiler CSV file ready for use with parse_db.py

For insertions and deletions we need to account for the strand when calculating
the position.

"""

import re
import difflib
import logging
from importlib import resources
import gffutils
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.Data.IUPACData import protein_letters
from typing import List


logger = logging.getLogger()


# defining regex patterns for each variant case
# note the WHO is using ! to symbolise stop codon rather than *
AA_CHANGE = re.compile(f'([{protein_letters}])([0-9]{{1,4}})([{protein_letters}!])')
ncRNA_CHANGE = re.compile(f'([acgt])([0-9]{{1,4}})([acgt])')
PROMOTER_CHANGE = re.compile(f"([acgt])(-[0-9]{{1,3}})([acgt])")
NUC_DELETION = re.compile(f"([0-9]{{1,4}})_(del)_([0-9]{{1,2}})_([actg]+)_([actg]+)")
NUC_INSERTION = re.compile(f"([0-9]{{1,4}})_(ins)_([0-9]{{1,3}})_([acgt]+)_([acgt]+)")
PROMOTER_DEL = re.compile(f"-([0-9]{{1,3}})_(del)_([0-9]{{1,4}})_([acgt]+)_([acgt]+)")
PROMOTER_INS = re.compile(f"-([0-9]{{1,3}})_(ins)_([0-9]{{1,4}})_([acgt]+)_([acgt]+)")


def get_gene_strands() -> dict:
    """
    Read the reference genome GFF file and return a dictionary with
    gene IDs as keys and the strand as the values
    Returns both CDS regions and rRNAs
    :return: dictionary
    """
    logger.info("Loading GFF file to get gene strands...")
    with resources.path('who2tbp.lib.data', 'genome.gff') as gff_file:
        db = gffutils.create_db(str(gff_file), ':memory:')
        return {rec.attributes['Name'][0]: rec.strand for rec in db.all_features()
                if rec.featuretype in ['gene', 'rRNA_gene'] and rec.attributes.get('Name', None) is not None}


def insertion_calc(ref_seq: str, alt_seq:str) -> tuple:
    """
    Given a reference sequence and an alternate sequence, identify
    the inserted bases and offset of the insertion relative to the
    initial base in the ref string
    :param ref_seq: the sequence in the reference genome
    :param alt_seq: the sequence in the mutated genome
    :return: offset_pos, inserted sequences
    """

    diff = difflib.ndiff(ref_seq.upper(), alt_seq.upper())
    insert_data = [(pos, base.strip("+").strip()) for pos, base in enumerate(diff) if base.startswith('+')]
    insert_seq = ''.join([rec[1] for rec in insert_data])
    offset = insert_data[0][0] - 1
    return offset, insert_seq


def promoter_offset(ref_seq: str, alt_seq: str) -> tuple:
    """
    Given a reference sequence to a promoter and an alternate sequence with a deletion,
    identify the position of the deletion in the reference sequence
    :param ref_seq: the reference sequence
    :param alt_seq: the alternate sequence
    :return: offset, del_len
    """
    diff = difflib.ndiff(ref_seq, alt_seq)
    del_data = [(pos, base.strip("-".strip())) for pos, base in enumerate(diff) if base.startswith("-")]
    offset = del_data[0][0] - 1
    len_del = len(del_data)
    return offset, len_del


def who2hgvs(var: str, this_gene_strands: dict) -> tuple:
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

    match_prom_del = PROMOTER_DEL.match(var)
    if match_prom_del:
        offset, del_len = promoter_offset(match_prom_del.group(4), match_prom_del.group(5))
        start_pos = int(match_prom_del.group(1)) + offset + 1
        end_pos = int(match_prom_del.group(3)) + start_pos - 1
        if end_pos - start_pos == 0:
            return gene, f"c.-{start_pos}del"
        else:
            return gene, f"c.-{start_pos}_-{end_pos}del"

    match_prom_ins = PROMOTER_INS.match(var)
    if match_prom_ins:
        indicator_pos = int(match_prom_ins.group(1))
        ref_seq = match_prom_ins.group(4)
        alt_seq = match_prom_ins.group(5)
        offset, ins_seq = insertion_calc(ref_seq, alt_seq)
        start_pos = indicator_pos + offset
        end_pos = start_pos + 1
        return gene, f"c.-{start_pos}_-{end_pos}ins{ins_seq}"


    # from here on end we MAY need strand information to make identify the
    # correct location of the changes
    strand = this_gene_strands.get(gene, None)
    if not strand:
        logger.warning(f"Did find strand information for {gene}. Assuming positive strand.")

    # deletions are given upstream from the indicated position
    # so, need to massage the output
    match_del = NUC_DELETION.match(var)
    if match_del:
        if strand == '+':
            start_pos = int(match_del.group(1)) + 1
            end_pos = int(match_del.group(1)) + int(match_del.group(3))
        else:
            start_pos = int(match_del.group(1)) - int(match_del.group(3))
            end_pos = int(match_del.group(1)) - 1
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
        offset, ins_seq = insertion_calc(ref_seq, alt_seq)
        if strand == '+':
            start_pos = indicator_pos + offset
            end_pos = start_pos + 1
        else:
            end_pos = indicator_pos - offset
            start_pos = end_pos - 1
        return gene, f"c.{start_pos}_{end_pos}ins{ins_seq}"


def who2tbdw(data: List[dict]) -> dict:
    """
    Given the output from WHO extraction, translate the variant code to 
    something the HGVS nomenclature

    :param data: 
    :return: 
    """
    gene_strands = get_gene_strands()

