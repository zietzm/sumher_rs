# `sumher_rs`

[LDAK](https://dougspeed.com/) is a powerful piece of software for genetics.
SumHer is a tool for estimating heritability and genetic correlation using GWAS summary statistics that is [implemented in LDAK](https://dougspeed.com/sumher/).

`sumher_rs` is a wrapper around LDAK that is optimized for computing many heritabilities and genetic correlations in parallel.
To do so, it avoids re-reading files from disk and avoids unnecessarily repeating checks and formatting.
This package is most suited to pre-aligned GWAS summary statistics.
It requires the same inputs as LDAK/SumHer, but can efficiently run for multiple phenotypes at once.

## Installation

If cargo is not installed, see [cargo installation](https://doc.rust-lang.org/cargo/getting-started/installation.html).

```bash
git clone --depth 1 https://github.com/zietzm/sumher_rs.git

cargo install --path sumher_rs
```

## Usage

`sumher_rs` is a command line tool.

To see a full list of commands and arguments, run

```bash
sumher_rs -h
```
