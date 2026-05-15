# `mantis-delta` preprint manuscript

This directory contains the source files for the `mantis-delta` research article
prepared for submission to **bioRxiv** as a preprint.

## Files

| File | Purpose |
|------|---------|
| `mantis_preprint.qmd` | Quarto source for the manuscript. Renders to PDF (XeLaTeX) or HTML. |
| `references.bib`     | BibTeX bibliography (Horn–Jackson, Feinberg, Goldbeter–Koshland, NUPACK, SciPy/NumPy/SymPy, etc.). |

## Rendering

```bash
# Install Quarto: https://quarto.org/docs/get-started/
cd paper
quarto render mantis_preprint.qmd --to pdf     # produces mantis_preprint.pdf
quarto render mantis_preprint.qmd --to html    # produces mantis_preprint.html
```

The PDF target uses XeLaTeX with the KOMA-Script `scrartcl` class, Times New
Roman, 11 pt, 1.25 line spacing, and 1-inch margins, matching the
typography requested by the bioRxiv submission system.

## Figures

The manuscript references PNG figures that already live in
`../examples/`:

- `cha_bifurcation.png`         — Fig. 2 (CHA dose-response)
- `cha_reaction_graph.png`      — Fig. 3 (CHA complex graph)
- `brusselator_oscillations.png` — Fig. 1 (Brusselator Hopf)
- `gk_switch_dose_response.png` — Fig. 4 (Goldbeter–Koshland switch)

To regenerate all figures from scratch:

```bash
python examples/04_cha_cascade.py
python examples/05_brusselator_chemostatted.py
python examples/06_goldbeter_koshland.py
```

## bioRxiv submission checklist

| Field (bioRxiv) | Source |
|---|---|
| Article category | YAML `article-category` |
| Subject area | YAML `subject-area` |
| Title | YAML `title` |
| External data URL | YAML `external-data` |
| Abstract | YAML `abstract` (also rendered into the PDF) |
| Author approvals | YAML `author-approvals` |
| Competing interests | YAML `competing-interests` |
| Author list / ORCID | YAML `author` |
| Funding | YAML `funding` |
| Distribution / reuse | YAML `distribution-reuse` |
| Manuscript file (PDF) | `quarto render ... --to pdf` |
| Image files (PNG) | `../examples/*.png` |
