# Decoder for hybrid toric-color (HTC) codes

Implementation of the decoding algorithm found in [arXiv:2011.xxxx](https://arxiv.org)

## Usage

- Use `gen_input.py` in input to generate a parameter file
  - `L` is linear lattice length for an LxL torus (must be a multiple of three)
  - `eP` is phase-flip error probability
  - `uP` is the morphing (unencoding) probability
  - `disorder` is the number of lattices to generate for each `uP`
- Use `monte_carlo.py` to run experiments using a parameter file
- Also included example script for running on a cluster (`array-job.sh`)

## Notes

To only morph red ball-like regions, set `rOnly` to `true` in `main.cpp`

## Software

This code makes use of [Blossom V](http://pub.ist.ac.at/~vnk/software.html)