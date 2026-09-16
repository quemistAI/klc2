#!/bin/bash
set -e
cd ~/klc2
mkdir -p raw/mn_atac
cut -f2 processed/mn_atac_selected.tsv | tail -n +2 | while read acc; do
  out="raw/mn_atac/${acc}.bigWig"
  if [ -s "$out" ]; then
    echo "have $acc"
  else
    echo "fetching $acc"
    curl -L -o "$out" "https://www.encodeproject.org/files/${acc}/@@download/${acc}.bigWig"
  fi
done
ls -lh raw/mn_atac/

