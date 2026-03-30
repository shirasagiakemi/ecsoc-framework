# Data Format & Acquisition

## Sample Data (included in this repository)

Two simulated RR interval files are provided for pipeline testing.
**These are NOT real patient data.**

| File | N | Description |
|------|---|-------------|
| `sample_rr_normal.csv` | 500 beats | Simulated healthy-like RR series (high variability, AR+noise) |
| `sample_rr_chf.csv`    | 500 beats | Simulated CHF-like RR series (low variability, high autocorrelation) |

### CSV Format

```
rr_ms,source
847.3,simulated_healthy
821.1,simulated_healthy
...
```

| Column | Unit | Description |
|--------|------|-------------|
| `rr_ms` | milliseconds | RR interval (beat-to-beat) |
| `source` | string | Data origin label |

---

## Real PhysioNet Datasets (not included — download separately)

All datasets used in the paper are openly available at [physionet.org](https://physionet.org/).

### Download via Python (wfdb)

```bash
pip install wfdb
```

```python
import wfdb
import numpy as np

# Example: NSR Healthy Database
record_list = wfdb.get_record_list('nsrdb')
for rec in record_list[:3]:
    ann = wfdb.rdann(f'nsrdb/{rec}', 'atr', pn_dir='nsrdb')
    rr_ms = np.diff(ann.sample) / ann.fs * 1000
    np.savetxt(f'data/raw/{rec}_rr.csv', rr_ms, header='rr_ms', comments='')
```

### Dataset Index

| PhysioNet ID | Cohort in paper | N | URL |
|---|---|---|---|
| `nsrdb` | NSR Healthy | 54 | https://physionet.org/content/nsrdb/ |
| `chf2db` | CHF2 (NYHA 1–3) | 29 | https://physionet.org/content/chf2db/ |
| `chfdb` | BIDMC CHF (NYHA 3–4) | 15 | https://physionet.org/content/chfdb/ |
| `cudb` | CUDB (acute VF) | 35 | https://physionet.org/content/cudb/ |
| `svtdb` | SVTDB (VT/VF) | 78 | https://physionet.org/content/svtdb/ |
| `mvedb` | MVEDB | 22 | https://physionet.org/content/mvedb/ |
| `cast-rr` | CAST | 1543 | https://physionet.org/content/cast-rr/ |

### Preprocessing Steps

1. Extract RR intervals from annotation files (`.atr`, `.qrs`)
2. Remove ectopic beats: `rr < 300 ms` or `rr > 2000 ms`
3. Remove artifacts: consecutive RR change > 20%
4. Minimum series length: **256 beats** (recommended: ≥ 512)
5. Save as CSV with `rr_ms` column header

```r
# R preprocessing example
preprocess_rr <- function(rr_raw) {
  rr <- rr_raw[rr_raw >= 300 & rr_raw <= 2000]          # range filter
  diff_pct <- abs(diff(rr)) / rr[-length(rr)]
  rr <- rr[c(TRUE, diff_pct <= 0.20)]                   # artifact filter
  rr
}
```

---

## Minimum Data Requirements for DFA

| Parameter | Minimum | Recommended |
|-----------|---------|-------------|
| Series length | 256 beats | ≥ 512 beats |
| Recording type | Any Holter | Long-term ≥ 20 h for Phase V assessment |
| Sample rate | Any (after RR extraction) | — |

---

## Privacy Note

**Never commit real patient data to this repository.**  
The `.gitignore` excludes `data/raw/`, `*.dat`, `*.hea`, and `*.atr` files.  
PhysioNet data are de-identified but should be handled per their respective licenses.
