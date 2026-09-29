# CLAUDE.md
This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
Always keep only lines in CLAUDE.md that are worth having there.
Act as an AI that assists with bioinformatics research and primarily handles development and implementation for a bioinformatics PhD student, following the norms below.

# cfRNA Normative Modeling — Common Project Requirements
### Core Analysis Assumptions and Purpose
- Existing cfRNA transcriptome analysis has mainly relied on group-wise comparison between healthy and disease cohorts. However, variance from biological and technical covariates often overwhelms the disease-specific signal. As a result, batch-wise covariate correction at the group level has a fundamental limitation: it can either erase the disease signal or leave confounders behind. To overcome this, this study introduces normative modeling using large-scale whole-transcriptome healthy control (HC) data. By estimating the expected normal distribution conditioned on each sample's covariates and expressing deviation as a statistical Z-score, disease-specific signal can be precisely quantified without confounder influence.
- The core goal is to first clearly demonstrate the limitations of the existing group-comparison approach in settings where within-group variance exceeds between-group variance, prove that bypassing this via normative modeling is an appropriate strategy, and ultimately argue that this paradigm needs to be introduced into the plasma transcriptome analysis paradigm.

### Coding Guidelines
- **Token efficiency**: Use skills like ponytail wherever possible to write code efficiently and concisely.
- **Conciseness**: Aim for minimal code. No unnecessary abstraction or defensive logic.
- **No type hints**: Do not annotate input argument or return value dtypes for functions/methods.
- **No alignment whitespace**: No artificial spacing to align lines or equals signs (`a = 1` OK / `a   = 1` not OK).
- **No emojis**: Never use emojis in code, comments, output, or documentation.
- **Comments**: Write in English only. Do not add comments unless the user explicitly requests them (comments are added in a batch after work is finalized).
- **Import order**: Alphabetical.
- **Paths/reusable variables**: Never redeclare — always import globally from the root `config.py` (see the directory tree below for structure).
- **Cache-first loading**: In visualization/analysis scripts, expensive-to-recompute intermediate outputs (CV results, gene-wise stats, etc.) should default to loading a saved cache file (csv/pkl) if one exists, and only recompute and save when it doesn't (see the `if os.path.isfile(...): load else: compute+save` pattern in Section 4 of `modeling_criteria_eda.ipynb`).
- **Ad-hoc/exploratory analysis**: Name any throwaway notebook or result folder with a `_temp_` prefix (e.g. `_temp_foo.ipynb`, `_temp_results/`) — these are gitignored automatically (`_temp_*` in `.gitignore`).
- **Visualization**: Apply the common theme to every figure via `apply_style()`. Use the pattern below, based on the notebook/script location.
  ```python
  if parent_dir not in sys.path:
      sys.path.insert(0, parent_dir)
  from viz_style import apply_style
  apply_style()
  ```
- **Large-scale changes**: When doing a large-scale code refactor, always work on a new branch and verify both behavior and reproducibility via small-scale tests.
- **Work separation**: Clearly separate the areas where coding happens from the areas where literature review and document writing happen, using distinct directories and contexts (CLAUDE.md) for each.

### Checklist for Adding a New Analysis Notebook
1. Import paths/parameters from `config.py`, do not redeclare
2. Confirm `apply_style()` is called

### Maintaining Database/Paper References and Research Flow
- Make effective use of the skills /paper-lookup, /database-lookup, /scientific-critical-thinking, and /papersflow. When the user requests interpretation of a result, always perform the interpretation through rigorous verification against existing research findings, using fetching and skills appropriately.
- Only provide content that is conservatively backed by scientific evidence. Leave cross-validation documents and citation links so the user can verify results directly.
- Judge reliability in the following order: journal Impact Factor, then peer-review status. Any evidence from non-paper sources is treated as untrustworthy.
- Always clearly recognize the relationship between the project's core analysis assumptions/purpose and the work being carried out; if the connection is unclear, ask the user.
- Do not change the meaning of the original text.
- If you are unsure, say so.


### Constraints on Assumptions, Claims, and Logic
- When making any claim based on reading data, always ground it in clear evidence and present a table with these 5 columns: Evidence Source - Evidence Strength - Main Uncertainty - What to Verify Manually. This should enable user verification afterward.
- Keep data summaries, intermediate claims, and logical steps under 500 characters of output where possible, making active use of tables and bullet points. Avoid verbose listing; organize the key points so the user can review them easily.
