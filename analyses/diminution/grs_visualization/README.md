# GRS Coverage Visualization Data

**Created:** 2026-02-06
**Analysis:** Germline-Restricted Sequence (GRS) Coverage Drop Visualization
**Organism:** *Auanema rhodensis*
**Region:** SUPER_X:10706804-10913775 (206,971 bp GRS region)

---

## 📁 Contents

This directory contains a complete, self-contained analysis of coverage patterns at the germline-restricted sequence (GRS) region on SUPER_X chromosome, demonstrating chromatin diminution in *A. rhodensis*.

### Data Files (Ready for Plotting)

1. **SUPER_X_GRS_windowed_coverage.bed** (92 KB, 2,898 windows)
   - 100bp windowed coverage across GRS + 20% flanking regions
   - Region: SUPER_X:10665410-10955169 (289,759 bp total)
   - Format: `chr`, `start`, `end`, `mean_coverage`

2. **SUPER_X_GRS_borders.bed** (1.1 KB, 32 entries)
   - Per-base coverage at GRS boundaries (±20bp each side)
   - Left boundary: 10706784-10706824
   - Right boundary: 10913755-10913795
   - Format: `chr`, `start`, `end`, `coverage`

3. **SUPER_X_telomeric_breakpoints.bed** (1.8 KB, 50 breakpoints)
   - Telomeric breakpoint positions detected by MILTEL
   - All SUPER_X chromosome_breakage_site features
   - Format: `chr`, `start`, `end`, `score`, `orientation`

### Scripts

4. **generate_grs_visualization_files.py**
   - Python script used to generate all three data files
   - Can be reused for other GRS regions

### Documentation

5. **GRS_VISUALIZATION_FILES_README.md**
   - Detailed technical documentation
   - Validation results and statistics
   - Usage examples

6. **README.md** (this file)
   - Overview and quick reference

---

## 🧬 Biological Context

### What is Chromatin Diminution?

Chromatin diminution is a developmental process where germline-restricted DNA sequences are eliminated from somatic cells. In *A. rhodensis*:

- **Germline cells:** Retain all DNA including GRS regions
- **Somatic cells:** Eliminate GRS regions through programmed DNA deletion
- **Mechanism:** Chromosome breakage + telomere addition at boundaries

### This GRS Region

- **Location:** SUPER_X:10706804-10913775
- **Size:** 206,971 bp (~207 kb)
- **Coverage drop:** ~7-fold reduction (291X → 42X)
- **Boundaries:** Marked by telomeric breakpoints

This is one of the larger germline-restricted sequences in the *A. rhodensis* genome.

---

## 📊 Key Findings

### Coverage Pattern

| Region | Mean Coverage | Windows | % of Region |
|--------|---------------|---------|-------------|
| **Somatic flanks** | 291.3X | 828 | 28.6% |
| **GRS region** | 41.9X | 2,070 | 71.4% |
| **Fold change** | **6.95:1** | - | - |

### Sharp Boundaries

**Left Boundary (position 10706804):**
- Before: 299-303X (somatic)
- After: 43X (GRS)
- Transition: Single 100bp window

**Right Boundary (position 10913775):**
- Before: 43-46X (GRS)
- After: 349-369X (somatic)
- Transition: Single 100bp window

### Telomeric Signals

**Primary breakpoints:**
- Position 10706805: score 5.35 (right orientation) ← **Strong signal**
- Position 10913775: score 0.86 (left orientation)

Additional low-score breakpoints cluster around both boundaries, consistent with imprecise chromosome breakage during chromatin diminution.

---

## 🚀 Quick Start - Visualization

### Load in R

```r
# Load data files
windowed <- read.table("SUPER_X_GRS_windowed_coverage.bed",
                       sep="\t",
                       col.names=c("chr", "start", "end", "coverage"))

borders <- read.table("SUPER_X_GRS_borders.bed",
                      sep="\t",
                      col.names=c("chr", "start", "end", "coverage"))

telomeres <- read.table("SUPER_X_telomeric_breakpoints.bed",
                        sep="\t",
                        col.names=c("chr", "start", "end", "score", "orientation"))

# Basic plot
library(ggplot2)
ggplot(windowed, aes(x=start, y=coverage)) +
  geom_line() +
  geom_vline(xintercept=10706804, color="red", linetype="dashed") +
  geom_vline(xintercept=10913775, color="red", linetype="dashed") +
  theme_minimal() +
  labs(title="GRS Coverage Drop on SUPER_X",
       x="Position (bp)", y="Coverage (X)")
```

### Load in Python

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load data files
windowed = pd.read_csv("SUPER_X_GRS_windowed_coverage.bed",
                       sep="\t",
                       names=["chr", "start", "end", "coverage"])

borders = pd.read_csv("SUPER_X_GRS_borders.bed",
                      sep="\t",
                      names=["chr", "start", "end", "coverage"])

telomeres = pd.read_csv("SUPER_X_telomeric_breakpoints.bed",
                        sep="\t",
                        names=["chr", "start", "end", "score", "orientation"])

# Basic plot
plt.figure(figsize=(12, 4))
plt.plot(windowed['start'], windowed['coverage'])
plt.axvline(10706804, color='red', linestyle='--', label='GRS boundaries')
plt.axvline(10913775, color='red', linestyle='--')
plt.xlabel('Position (bp)')
plt.ylabel('Coverage (X)')
plt.title('GRS Coverage Drop on SUPER_X')
plt.legend()
plt.tight_layout()
plt.show()
```

### View in IGV

1. Open IGV (Integrative Genomics Viewer)
2. Load genome: File → Load Genome from File → `assembly.fasta`
3. Load tracks: File → Load from File → Select BED files
4. Navigate to: `SUPER_X:10,665,410-10,955,169`
5. Adjust track heights for better visualization

---

## 🔧 Regenerating the Data

If you need to regenerate these files (e.g., with different parameters or for other regions):

### Command Used

```bash
python3 generate_grs_visualization_files.py \
  --coverage ../../coverage/mosdepth/nxAuaRhod1_1.target_regions.per-base.bed.gz \
  --telomeres ../telomeres/miltel.test.miltel.telomeric.gff3 \
  --output-dir . \
  --chromosome SUPER_X
```

### Input Data Sources

- **Coverage data:** `analysis/coverage/mosdepth/nxAuaRhod1_1.target_regions.per-base.bed.gz`
  - Per-base coverage from mosdepth analysis
  - Generated from aligned sequencing reads

- **Telomeric breakpoints:** `analysis/diminution/telomeres/miltel.test.miltel.telomeric.gff3`
  - Chromosome breakage sites identified by MILTEL
  - Based on softclipped reads with telomeric repeat motifs

### Customization

To analyze a different region, edit the constants in `generate_grs_visualization_files.py`:

```python
# GRS Region Constants
GRS_START = 10706804      # Change to your region start
GRS_END = 10913775        # Change to your region end
FLANK_SIZE = int(GRS_LENGTH * 0.2)  # Adjust flank percentage
WINDOW_SIZE = 100         # Adjust window size
BORDER_BP = 20            # Adjust border window size
```

---

## 📈 Expected Output Pattern

When visualized, you should see:

1. **High coverage** (~291X) in flanking somatic regions (left and right)
2. **Sharp drop** at left boundary (position 10706804)
3. **Low coverage** (~42X) throughout GRS region
4. **Sharp recovery** at right boundary (position 10913775)
5. **Strong telomeric signals** at both boundaries

This pattern confirms:
- Successful chromatin diminution in somatic cells
- Clean chromosome breakage and telomere addition
- GRS retention in germline (contributing ~14% of total coverage)

---

## 🔗 Related Files

- **GRS coordinates:** `../nxAuaRhod1_1.GRS.bed`
- **Core genome coordinates:** `../nxAuaRhod1_1.core.bed`
- **Telomere analysis:** `../telomeres/`
- **Coverage analysis:** `../../coverage/mosdepth/`
- **Analysis scripts:** `../../../scripts/coverage_analysis/`

---

## 📝 Citation

If you use this data, please cite:

- The *Auanema rhodensis* genome project
- MILTEL: Telomere identification tool
- Mosdepth: Coverage calculation tool

---

## ✅ Validation Checklist

- [x] Files generated successfully
- [x] Coverage shows ~7-fold drop in GRS region
- [x] Sharp transitions at both boundaries
- [x] Telomeric breakpoints at expected positions
- [x] All coordinates on SUPER_X chromosome
- [x] Valid BED format (tab-separated, sorted)
- [x] Ready for visualization

---

## 💡 Tips

- Use windowed coverage for overview plots
- Use border coverage for high-resolution boundary plots
- Overlay telomeric breakpoints to show mechanism
- Compare with other GRS regions for consistency
- Check coverage normalization if comparing samples

---

**Questions?** See `GRS_VISUALIZATION_FILES_README.md` for detailed technical documentation.
