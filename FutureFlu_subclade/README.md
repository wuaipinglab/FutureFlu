# FutureFlu_subclade

[简体中文](README.zh.md)

Extension of FutureFlu for A/H1N1pdm09, A/H3N2, and B/Victoria in seasons
2025–2026 (north and south). 2025 uses clade labels; only 2026 is evaluated
at subclade resolution.

## Layout

- `scripts/` — analysis entry points
- `data/` — subclade definitions, counts, HA1 clade priors, configs, positivity, EVEscape
- `outputs/predictions/` — prediction tables
- `raw_inputs/` — place private inputs here for a full rerun

## Run

```bash
pip install -r ../requirements.txt
python scripts/run_pipeline.py definitions
python scripts/run_pipeline.py counts
python scripts/run_pipeline.py aux
```

`definitions` / `counts` refresh subclade labels and counts.
`aux` scores subclade risk components against truth and builds E/G/D combination
tables under `outputs/predictions/`.
