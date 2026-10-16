# Test fixtures

`data-raw/trim_fixtures.R` builds these files. Each `.spf.xz` file holds the
bytes that `read_spf()` uses for the selected reads of a published Sequencher
project, with every other byte set to zero, compressed with xz. The trimmed
files are not valid Sequencher projects. `expected.rds` holds the values that
`read_spf()` returns for the same reads from the untrimmed projects.

| Fixture | Source | Reads |
|---|---|---|
| `S1-trimmed.spf.xz` | S1 Sequencher File | `5C-est2-6-M13-FP.ab1` |
| `S8-trimmed.spf.xz` | S8 Sequencher File | `MC-1A-2-M13F(-47).ab1`, `MC-1A-37-M13F(-47)_R.ab1` |

Both source files are supplementary data to Claussin C, Chang M (2016).
Multiple Rad52-mediated homology-directed repair mechanisms are required to
prevent telomere attrition-induced senescence in *Saccharomyces cerevisiae*.
*PLOS Genetics* 12(7): e1006176. <https://doi.org/10.1371/journal.pgen.1006176>,
distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The trimmed versions are modified from the originals.

The test of all published projects also uses the six projects in
[Zenodo record 5557681](https://zenodo.org/records/5557681) (CC0). It runs
only when `SEQSALVAGE_TEST_DATA` names a directory holding some of the 14
projects, which you download yourself from the URLs in `helper-published.R`.
Each file is checked against its published md5 checksum.
