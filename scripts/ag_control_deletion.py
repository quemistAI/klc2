#!/usr/bin/env python3
"""Control: a matched 217-bp deletion elsewhere in the same lncRNA intron.
If the lncRNA score is similar, the effect is an artifact of deleting within
the gene body, not regulation by the SPOAN element."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env
API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# same intron, ~2 kb upstream of the element, outside the promoter/CpG shore.
# CHECK THIS INTERVAL IN UCSC FIRST: still intronic, no cCRE, no CpG island.
CTRL_START, CTRL_END = 66254800, 66255016      # 217 bp

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

import requests  # not in this env; see note below
