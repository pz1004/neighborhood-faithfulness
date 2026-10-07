#!/usr/bin/env bash
# Stage 2 -- the recordings, 6.5 GB: the preprocessed THINGS-EEG2 test files of 10
# participants, and the 10 Allen two-photon experiments that docs/note_protocol_allen.md,
# Section 2, selects. Each fetcher records the SHA-256 digest of every file in
# data/*/manifest.json; this rewrites those tracked manifests and the fetch logs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 analyses/fetch/fetch_things_eeg2_test.py > data/fetch.log
curl -f -o data/allen_2p/ophys_experiments.json \
  https://allen-brain-observatory.s3.us-west-2.amazonaws.com/visual-coding-2p/ophys_experiments.json
python3 analyses/fetch/fetch_allen_2p.py > data/allen_2p/fetch.log
