# who2tbp: Convert WHO TB antimicrobial resistant variants to TBProfiler format

## Installation

```bash
pip install git+https://github.com/MDU-PHL/who2tbp.git
```

## Use

Get all the high confidence variants associated with resistance:
```bash
who2tbp WHO-UCN-GTB-PCI-2021.7-eng.xlsx > tbdb_who.py
```

Get all the variants only interim associated with resistance:
```bash
who2tbp --filter assoc_resistance_interim WHO-UCN-GTB-PCI-2021.7-eng.xlsx > tbdb_who_assoc_interim.py
```

Get help:
```bash
who2tbp -h 
```

## Complete usage

```bash
usage: Convert WHO Excel sheet with MTB mutations to TBProfiler database format [-h] [-f {assoc_resistance,no_assoc,assoc_resistance_interim,no_assoc_interim,combo,uncert_signif,all}] [-o OUTFILE] INFILE

positional arguments:
  INFILE                The WHO Excel sheet

optional arguments:
  -h, --help            show this help message and exit
  -f {assoc_resistance,no_assoc,assoc_resistance_interim,no_assoc_interim,combo,uncert_signif,all}, --filter {assoc_resistance,no_assoc,assoc_resistance_interim,no_assoc_interim,combo,uncert_signif,all}
                        Limit to single category (default: assoc_resistance)
  -o OUTFILE, --outfile OUTFILE
```


## Develop

1. Fork the repository:

```bash
gh fork repo MDU-PHL/who2tbp
```

2. Install all dependencies:

```bash
cd who2tbp
conda create -n who2tbp -f requirements.txt
```

3. Make changes

4. Commit changes to your fork

5. Create a pull request

Check out this guide on creating pull requests if you are unsure how to: https://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR

