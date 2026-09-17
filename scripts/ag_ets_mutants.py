#!/usr/bin/env python3
"""Point-mutate the two ETS cores in the corrected SPOAN element and score
predicted effect on KLC2. Deletions change spacing; substitutions do not,
so this isolates the binding site from the geometry.
RUN IN: conda env `alphagenome`."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env
API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

FA, PAD = "raw/element_TRUE_flank500_hg38.fa", 500
EL_START = 66257086                  # hg38, 1-based

CORE1 = (109, 117)                   # element coords, 1-based inclusive
CORE2 = (168, 176)
WT_CORE = "ACCGGATGT"
# GGA -> TTA destroys the ETS GGAA/T core while preserving length and GC
MUT_CORE = "ACCTTATGT"

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()

for label, (a, b) in (("core1", CORE1), ("core2", CORE2)):
    obs = seq[PAD + a - 1: PAD + b]
    assert obs == WT_CORE, f"{label}: expected {WT_CORE}, found {obs}"
    print(f"{label} verified at element {a}-{b}: {obs}")

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]


def substitution(a, b, new):
    """Variant replacing element positions a..b with `new` (same length)."""
    g = EL_START + a - 1
    ref = seq[PAD + a - 1: PAD + b]
    assert len(ref) == len(new)
    return genome.Variant(chromosome="chr11", position=g,
                          reference_bases=ref, alternate_bases=new)


def deletion(a, b):
    g = EL_START + a - 1
    return genome.Variant(chromosome="chr11", position=g - 1,
                          reference_bases=seq[PAD + a - 2] + seq[PAD + a - 1: PAD + b],
                          alternate_bases=seq[PAD + a - 2])


VARIANTS = {
    "full_217_del":  deletion(1, 217),
    "core1_mut":     substitution(*CORE1, MUT_CORE),
    "core2_mut":     substitution(*CORE2, MUT_CORE),
    "core1_del":     deletion(*CORE1),
    "core2_del":     deletion(*CORE2),
}

frames = []
for label, v in VARIANTS.items():
    print(f"scoring {label} ...")
    s = model.score_variant(
        interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=v, variant_scorers=[scorer])
    df = variant_scorers.tidy_scores(s)
    df["edit"] = label
    frames.append(df)

# both cores mutated at once: two separate variants can't be combined in one
# score_variant call, so this is noted as a limitation rather than run here.
print("\nNOTE: the double mutant requires predict_sequence on a custom sequence,"
      "\nnot score_variant. Deferred - see notebook.")

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_ets_mutant_scores.tsv", sep="\t", index=False)

FIB = r"fibroblast|BJ|IMR-90|HFFc6|AG04450"
NERVE = r"tibial nerve|sciatic nerve|spinal cord|motor neuron"
pd.set_option("display.width", 200)
for name, pat in (("fibroblast", FIB), ("nerve/cord", NERVE)):
    sub = out[out.biosample_name.astype(str).str.contains(pat, case=False, na=False)]
    piv = (sub[sub.gene_name.isin(["KLC2", "ENSG00000255320"])]
           .groupby(["edit", "gene_name"]).raw_score.mean().unstack().round(4))
    print(f"\n=== {name} mean score ===")
    print(piv.to_string())
