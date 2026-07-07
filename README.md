# Separating spread from detection — effort-corrected freshwater invasion dynamics from GBIF

A complete, reproducible study — **code + manuscript** — built to be submittable to a
Q1 journal (*Diversity and Distributions* / *NeoBiota* / *Biological Invasions*), using
**only GBIF data**.

## The idea in one paragraph

"GBIF-only" rules out standard SDM (no WorldClim/CHELSA predictors), so the study is
designed around what GBIF alone can answer rigorously: the **temporal dynamics** of
invasions. The novel, timely angle is that raw occurrence accumulation conflates two
things — genuine range expansion and the recent explosion of sampling effort (the
citizen-science boom). We quantify that bias across 14 European freshwater invaders and
correct it using **three independent methods** built from a sampling-effort denominator
taken from the same database: continental effort rarefaction (target-group background),
the same rarefaction restricted to each species' accessible area (Barve et al. 2011), and
a list-length detection model (Isaac et al. 2014). Comparing them shows the correction is
not an artefact of any one method's assumptions — a robustness check most GBIF-based
spread-rate papers skip entirely, which (on top of the effort-correction gap itself) is
what makes it Q1-worthy.

## ⚠️ Honesty note — read first

This environment **cannot download from GBIF** (network is restricted to a few package
mirrors), so I could not run the analysis and produce real numbers here. What you have is:

- a **fully runnable pipeline** (`R/`) — run it with your GBIF login and it produces every
  figure, table, and the numbers;
- a **complete manuscript** (`manuscript/manuscript.md`) — Introduction, Methods, and
  Discussion are finished prose; the **Results section is a template** with `⟨tokens⟩`
  that populate from the pipeline output.

I deliberately did **not** invent results and dress them as findings — that would be
fabrication. Run the pipeline, then the `⟨tokens⟩` become real. Nothing in the manuscript
should be reported until it is backed by an actual run.

## File map

```
gbif-invasion-tempo/
├── README.md                    ← this file
├── R/
│   ├── 00_setup.R               config, packages+versions, species table, helpers
│   ├── 01_data_acquisition.R    GBIF downloads (focal + target group) + cleaning   [internet]
│   ├── 02_range_dynamics.R      gridding, naive + 3 correction methods, phases, spread rate,
│   │                            cross-method agreement (Spearman rho, Kendall's W)
│   ├── 03_synthesis_figures_tables.R   Figs 1–6, Tables 1–3, results_values.json
│   └── run_all.R                sources everything in order
├── manuscript/
│   └── manuscript.md            full draft; Results is a token template
├── output/
│   ├── figures/                 Fig1–Fig6 (300 dpi, EPSG:3035, colourblind-safe)
│   ├── tables/                  Table1–3 + spread_rates / bias_summary / method_agreement / DOIs / etc.
│   └── results_values.json      numbers injected into the manuscript
└── data/                        cached GBIF downloads + cleaned .rds
```

## Correction methods compared

The pipeline reconstructs range dynamics under a naive (uncorrected) model and three
independent effort corrections, so the manuscript can report whether "naive reconstructions
are biased" is robust to the choice of correction rather than an artefact of one method:

| Method                 | Paradigm                        | Effort denominator          | Where |
|-------------------------|---------------------------------|------------------------------|-------|
| `naive`                 | none (baseline being tested)    | —                             | `naive_aoo()` |
| `rarefied`              | resampling to constant effort   | continental target group      | `rarefied_aoo()` |
| `rarefied_regional`     | resampling to constant effort   | accessible area (150 km buffer, Barve et al. 2011) | `restrict_target_accessible()` + `rarefied_aoo()` |
| `listlength`            | detection-probability modelling | list-length GLM (Isaac et al. 2014; Szabo et al. 2010) | `listlength_aoo()` |

`accessible_buffer_km` in `R/00_setup.R` (`cfg`) controls the accessible-area radius.

## How to run

1. Install R 4.5.0 and the packages listed at the top of `R/00_setup.R`.
2. Put GBIF credentials in `~/.Renviron` and restart R:
   ```
   GBIF_USER=yourname
   GBIF_PWD=yourpassword
   GBIF_EMAIL=you@institution.org
   ```
3. From the project root:
   ```bash
   Rscript R/run_all.R
   ```
   Step 01 is the slow part — GBIF prepares each download asynchronously and the
   target-group sets (e.g. all European Actinopterygii) are large. Downloads are cached,
   so re-runs only repeat the analysis (minutes).
4. Fill the manuscript: replace each `⟨token⟩` with the matching value from
   `output/results_values.json` and `output/tables/`.

## Token → source mapping (Results section)

| Token(s) in manuscript                              | Source |
|-----------------------------------------------------|--------|
| `⟨n_species⟩`, `⟨n_records_total⟩`, `⟨n_records_median⟩`, `⟨year_range⟩`, `⟨grid_km⟩` | `output/results_values.json` |
| `⟨median_spread_naive⟩`, `⟨median_spread_rarefied⟩` | `results_values.json` / `Table3_spread.csv` |
| `⟨median_spread_regional⟩`, `⟨median_spread_listlength⟩` | `results_values.json` / `Table3_spread.csv` |
| `⟨median_bias_ratio⟩`                                | `results_values.json` / `bias_summary.csv` |
| `⟨bias_cor_r⟩`, `⟨bias_cor_ci⟩`, `⟨bias_cor_p⟩`     | `results_values.json` |
| `⟨kendall_w⟩`, `⟨agreement_rho_min⟩`, `⟨agreement_rho_max⟩` | `results_values.json` / `method_agreement.csv` |
| `⟨jackknife_r_min⟩`, `⟨jackknife_r_max⟩`, `⟨jackknife_p_min⟩` | `results_values.json` / `jackknife_bias_cor.csv` |
| `⟨jackknife_r_drop_dbugensis⟩`, `⟨jackknife_r_drop_amelas⟩` | `results_values.json` / `jackknife_bias_cor.csv` |
| `⟨fastest_species⟩`                                  | `results_values.json` |
| lag-end / saturation-onset values                    | `Table2_phases.csv` |
| per-species spread + CIs                             | `Table3_spread.csv` / Fig. 4 |
| geographic / temporal descriptions                   | read off Fig. 1 and Fig. 2 |

## Easy swap points (make it your study)

- **Species / groups / target-group proxy** → `species_tbl` in `R/00_setup.R`.
- **Region** → `cfg$bbox` (and map limits in Fig. 1); switch CRS if you leave Europe
  (EPSG:3035 is Europe-specific; use another equal-area CRS elsewhere).
- **Grid resolution / temporal window** → `cfg$cell_km`, `cfg$year_min/max`.
- **Journal** → drafted for *Diversity and Distributions*; the structure fits *NeoBiota*
  and *Biological Invasions* with minor formatting changes.

## Design choices baked in (per your profile)

- Equal-area analysis (EPSG:3035, ETRS89-LAEA); area-based metrics require it.
- `set.seed()`, pinned package versions, explicit namespacing, native pipe.
- Uncertainty on every estimate (bootstrap percentile intervals; regression CIs).
- Colourblind-safe palettes (Okabe–Ito / viridis); 300 dpi; scale bars on maps.
- Effort correction compared across three independent methods (continental rarefaction,
  accessible-area rarefaction, list-length model), with dynamic occupancy models and the
  reference Frescalo (`sparta`) implementation flagged as heavier alternatives beyond
  those three.

## Sensible next steps

- Run it, then send me `results_values.json` + the figures and I'll (a) fill the Results
  prose, (b) tighten the Discussion validation paragraphs against the real numbers, and
  (c) convert the manuscript to **.docx** or a **parameterized .Rmd** that auto-populates.
- Optional rigor upgrades: dynamic occupancy models (`unmarked`) and the reference Frescalo
  implementation (`sparta`, benchmarked against a full community matrix rather than a
  single target-group ratio) remain heavier alternatives beyond the three corrections and
  the leave-one-species-out jackknife already implemented.
```
