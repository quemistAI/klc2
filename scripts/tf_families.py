"""Shared TF family classification. Imported by other scripts so the element
and the background are classified by identical rules."""

import re

FAMILIES = [
    ("ETS",      r"(ETS|ETV|ELK|ERF|FEV|FLI|ERG|GABPA|EHF|ELF[0-9]|SPI1|SPDEF)"),
    ("ZNF_C2H2", r"(ZNF|Zfp|ZKSCAN|ZSCAN|ZBED|PATZ|WT1|FEZF|BCL11|IKZF|ZIC|EGR|KLF|SP[0-9])"),
    ("ZBTB",     r"(ZBTB|Zbtb)"),
    ("NFKB",     r"(NFKB|REL[AB]?$)"),
    ("INSM",     r"INSM"),
    ("CTCF",     r"CTCF"),
    ("NR",       r"(NR[0-9]|RXR|RAR|THRB|ESR|PPAR|VDR)"),
    ("HOMEO",    r"(HOX|MEIS|PBX|TLX|DRGX|FOXI)"),
]

# PRESPECIFIED 2026-09-17, before running the background comparison:
# families counted as repressor-associated for the enrichment test.
# C2H2/KRAB zinc fingers and BTB-POZ factors are predominantly repressors;
# INSM1 is an established transcriptional repressor. NFKB excluded (p50/p52
# homodimers repress but heterodimers activate - context dependent). ETS
# excluded (predominantly activators, though ERF and ETV3 are repressors).
REPRESSOR_FAMILIES = {"ZNF_C2H2", "ZBTB", "INSM"}


def family(motif_name):
    for name, pat in FAMILIES:
        if re.search(pat, str(motif_name), re.I):
            return name
    return "other"
