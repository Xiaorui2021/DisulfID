# DisulfID

**Computer-assisted interpretation tool for disulfide-connectivity assignment in double-loop cysteine-rich peptides by enzymatic digestion and LC-MS/MS.**

DisulfID enumerates all tryptic-digestion combinations across the three bridge hypotheses of a four-cysteine peptide (**beads**, **globular**, **ribbon**), builds a theoretical b/y fragment library, and for each MS2 spectrum proposes connectivity candidates supported by **loop ions** (carrying two cysteines within a single chain) and **single-cysteine b-/y-linking ions** (indicating cross-chain pairing).

## Bridge hypotheses

| Name      | Pairing      | Cysteine arrangement |
|-----------|--------------|----------------------|
| beads     | C1–C2 / C3–C4 | adjacent            |
| globular  | C1–C3 / C2–C4 | crossover           |
| ribbon    | C1–C4 / C2–C3 | nested              |

## Features

- Theoretical fragment library enumerating all digestion combinations × three bridge hypotheses
- Both **b and y ions** participate in topology discrimination (loop ions + single-cysteine linking ions)
- MGF parsing with precursor isotope offset compensation (n = 0/1/2)
- Weighted scoring: `score = 100 × (0.7 × weighted ion coverage + 0.3 × intensity coverage)`
- Peak-dedup view (one representative spectrum per RT/mass cluster) and all-scan view
- Annotated MS2 figures (PNG + vector PDF)
- CSV outputs (`disulfide_judgment.csv`, `disulfide_judgment_allscans.csv`)
- HTML and Markdown reports
- PyQt6 graphical user interface; command-line mode also available

## Installation

### Option 1: Windows executable

Download `disulfide_analyzer.exe` from the [Releases](../../releases) page and run directly — no Python installation required.

### Option 2: Run from source

```bash
git clone https://github.com/[your-org]/DisulfID.git
cd DisulfID
pip install -r requirements.txt
python disulfide_analyzer.py
```

## Requirements

- Python 3.10+
- PyQt6
- matplotlib
- numpy

See `requirements.txt` for pinned versions.

## Usage (GUI)

1. Launch `disulfide_analyzer.exe` (or `python disulfide_analyzer.py`).
2. Enter the peptide sequence (four cysteines) and advanced parameters if needed.
3. Select the MGF file exported from your LC-MS/MS software.
4. Click **Run** — results appear in the table; switch between peak-dedup and all-scan views.
5. Click any row to inspect the annotated MS2 spectrum and evidence list.
6. Export CSV / HTML / Markdown reports via the **Export** button.

## Usage (command line)

```bash
python disulfide_analyzer.py input.mgf -o output_dir --sequence ACQEFCTYERDCASFDFDPSSNKC
```

## Scope and limitations

- Designed for **single-chain peptides with exactly two disulfide bonds** (four cysteines).
- The tool proposes per-scan candidates; the final intact-peptide connectivity should be assigned by the analyst from the consensus of fragment-level calls cross-checked against the major MS1 species.
- Spectra dominated by a short terminal fragment may yield degenerate (globular/ribbon indistinguishable) calls — the software flags these explicitly.
- Extension to three or more disulfide bonds is not supported in the current version.

## Citation

If you use DisulfID in your research, please cite:

> Bi Q, et al. *DisulfID: Computer-Assisted Disulfide Connectivity Assignment in Two-Disulfide-Bridge Peptides by Enzymatic Digestion and LC-MS/MS.* Analytical Chemistry, 2026.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contact

- Qirui Bi — bqrui2010@126.com
Shanghai Academy of International Standardization for Traditional Chinese Medicine, Shanghai Research Center for Modernization of Traditional Chinese Medicine, Shanghai University of Traditional Chinese Medicine, Shanghai 201203, China
