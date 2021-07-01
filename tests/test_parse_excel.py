"""
Test parse excel
"""

import pytest

from who2tbp.lib.parse_excel_data import parse_file


@pytest.mark.parametrize("conf_class,exp_len", [('assoc_resistance', 196),
                                                ('assoc_resistance_interim', 1004),
                                                ('uncert_signif', 15910),
                                                ('no_assoc_interim', 33),
                                                ('no_assoc', 213),
                                                ('combo', 40)])
def test_parse_file(conf_class, exp_len):
    """
    Read the file and some filters, and return a dictionary of the association
    :return:
    """
    with open("tests/data/WHO-UCN-GTB-PCI-2021.7-eng.xlsx", 'rb') as filehandle:
        data = parse_file(filehandle, this_filter=conf_class)
    assert len(data) == exp_len
