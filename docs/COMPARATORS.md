# Comparator sources for the HemaClone paper

Existing clonal haematopoiesis (CH) calling tools and workflows to compare
against. Details marked *to verify* come from abstracts or summaries and must
be checked against the full papers before they enter the manuscript.

---

## 1. ArCH

**Citation.** Chan ICC, Panchot A, Schmidt E, McNulty S, Wiley BJ, Liu J,
Turner K, Moukarzel L, Wong WSW, Tran D, Beeler JS, Batchi-Bouyou AL, Machiela
MJ, Karyadi DM, Krajacich BJ, Zhao J, Kruglyak S, Lajoie B, Levy S, Patel M,
Kantoff PW, Mason CE, Link DC, Druley TE, Stopsack KH, Bolton KL. ArCH:
improving the performance of clonal hematopoiesis variant calling and
interpretation. *Bioinformatics*. 2024;40(4):btae121.
doi:10.1093/bioinformatics/btae121

**Links.** Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11014783/ ·
Code: https://github.com/kbolton-lab/ArCH

**What it does**
- End-to-end, cloud-based workflow for CH SNVs and short indels.
- Combines four variant callers (which ones: *to verify*) and filters calls on variant characteristics and sequencing error-rate estimates.
- Customizable parameters for different sequencing technologies and datasets.
- Main target: deep targeted sequencing.

**Validation.** Six AML tumor:normal dilution series, 31 blood samples with
orthogonal validation, 26 blood samples with technical replicates. Reported
improved sensitivity and PPV for low-VAF CH variants over standard calling
(numbers: *to verify*).

**Relevance for the comparison.** The most direct comparator for calling
performance. Its ensemble calling and error-rate model target exactly the
low-VAF artefact problem where HemaClone is weakest (recurrent indels in
low-complexity regions). Whether it runs without a matched normal and on WES:
*to verify*.

---

## 2. qcCHIP

**Citation.** Liu X, Tang YH, Blachly J, Edge S, Jakubek YA, McCarter M, Naqash
AR, Nepple KG, Osman A, Reilley MJ, Riedlinger G, Salhia B, Schneider BP,
Shriver C, Churchman ML, Rounbehler RJ, Teer JK, Gillis N, Teng M. qcCHIP: an R
package to identify clonal hematopoiesis variants using cohort-specific data
characteristics. *Bioinformatics*. 2025;41(9):btaf522.
doi:10.1093/bioinformatics/btaf522. PMID: 40973027

**Links.** Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12466928/ ·
Code: https://github.com/tenglab/qcCHIP · Archive: doi:10.5281/zenodo.16421861

**What it does**
- Post-call filtering only: no alignment or calling. Input is single-sample VCFs (Mutect2 given as example), recommended ANNOVAR-annotated.
- Functions: `vcf2input` (merge and simplify VCFs), `CHIPfilter` (apply filters), `qcCHIP` (permutation analysis to guide threshold choice).
- Filter classes:
  - technical: VAF, DP, strand odds ratio, SAF/SAR, TLOD
  - functional: effect, repetitive-region location, frequency in gnomAD, dbGaP, COSMIC
  - individual: comparison with a paired non-blood sample, optional
  - population: allele- and locus-level prevalence in the cohort
- Modes: blood only, or blood paired with a non-blood control.
- Data: WES of blood (with or without paired tumor) and ultra-deep (~5000x) error-corrected targeted sequencing.

**Validation.** ORIEN breast (n=1338), lung (n=567), colorectal (n=1132), and an
ultra-deep breast cohort (n=406); permutation consistency compared with
manually curated CH calls (precision and recall). Reported cutoffs: VAF about
0.02–0.025 for WES and about 0.01 for ultra-deep data; prevalence 0.1; DP
about 20; SOR about 2–3; SAF/SAR about 3–5.

**Limitations noted.**
- Defaults tuned for typical WES.
- Prevalence consistency unstable with a few high-prevalence, low-VAF mutations.
- SAF/SAR did not converge in the ultra-deep cohort.
- `qcCHIP` evaluates only five metrics.
- The underlying data are available on request only.

**Relevance for the comparison.** Complementary rather than competing: it
filters calls that HemaClone produces. Its strand-bias and repetitive-region
filters are the ones HemaClone lacks. Its population (within-cohort prevalence)
filter has no counterpart by design: HemaClone classifies each sample from its
own data only. A natural
analysis is HemaClone calls with and without qcCHIP filtering.

---

## 3. CloneSifter (not a comparator)

**Citation.** Feldman D, Tsai F, Garrity AJ, O'Rourke R, Brenan L, Ho P,
Gonzalez E, Konermann S, Johannessen CM, Beroukhim R, Bandopadhayay P, Blainey
PC. CloneSifter: enrichment of rare clones from heterogeneous cell populations.
BioMed Central, 2021. https://hdl.handle.net/1721.1/131939 (journal and DOI:
*to verify*)

**Why it is listed.** Only as the name collision that ruled out "CloneSift"
during naming. It is a wet-lab method for isolating rare cell clones, not a CH
caller, and does not belong in the performance comparison.

---

## Suggested additional comparator

**Vlasschaert C, et al.** A practical approach to curate clonal hematopoiesis
of indeterminate potential in human genetic data sets. *Blood*.
2023;141(18):2214–2223.
- Record: https://stacks.cdc.gov/view/cdc/231491
- Stepwise filtering on sequencing metrics, variant annotation and population-based associations, with gene-specific rules (e.g. TET2, ASXL1).
- Applied to UK Biobank exomes and All of Us genomes (about 550,000 individuals).
- The de facto reference for tumor-only CHIP curation in large WES cohorts. Caller, thresholds, whitelist and code availability: *to verify* in the full paper.

## Positioning of HemaClone (draft)

| | HemaClone | ArCH | qcCHIP | Vlasschaert 2023 |
|---|---|---|---|---|
| Scope | FASTQ to CH status | end-to-end | post-call filter | curation protocol |
| Callers | Mutect2 | 4-caller ensemble | any (VCF input) | *to verify* |
| Matched normal needed | no | *to verify* | optional | no |
| Main data type | WES | deep targeted | WES, ultra-deep targeted | WES, WGS |
| Artefact model | Mutect2 filters, orientation, contamination, external PoN | error-rate estimates | strand bias, repeats, cohort prevalence, permutation | sequencing metrics, gene-specific rules |
| Deployment | local Docker, hard CPU and RAM limits, batch control | cloud | R package | protocol |
| Published benchmark | none yet | yes | yes | yes |

Claimed strengths to substantiate: end-to-end reproducibility from FASTQ,
resource-bounded local execution, per-criterion auditable outputs,
maftools-ready MAFs including CH negatives, preflight and fail-loud checks.
Gaps to close or acknowledge: single caller, no strand-bias or repeat-region
filters, no benchmark. Design difference to explain: no within-cohort
recurrence filter, because each sample is classified independently.
