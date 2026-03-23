# Kidney & Tumor Segmentation Post-processing

Post-processing tool for kidney + tumor + cyst segmentation labels.

> 한국어 버전: [README.md](README.md)

## Folder Structure

```
segmentation/
├── S003/
│   ├── S003_Segmentation_A.nii.gz    # Segmentation (A phase)
│   ├── S003_image_A.nii.gz           # CT image (A phase)
│   └── backup_original/              # Original backup (auto-created)
├── S004/
│   ├── S004_segmentation_A.nii.gz
│   ├── S004_image_A.nii.gz
│   └── ...
└── segtools.py
```

### Label Definitions

| Label | Name | Description |
|-------|------|-------------|
| `0` | Background | Outside organ |
| `1` | Kidney | Kidney parenchyma |
| `2` | Tumor | Tumor within kidney |
| `3` | Cyst | Cyst within kidney |

---

## Usage

```bash
python segtools.py <case_folder>
```

```bash
# Example
python segtools.py S004
```

### Workflow

```
1. Select phase (A / D / P / all)
2. Select function number
3. Enter sub-choices and parameters (press Enter for defaults)
4. Execute → auto backup → auto save
5. Review result, then enter (continue) or r (rollback)
6. Proceed to next function or q to quit
```

- Enter `q` at any input step to cancel the current function.
- Enter `r` at any input step to rollback to the previous state.
- On first run, the original file is automatically backed up to `backup_original/`.

---

## Function List

| # | Function | Description |
|---|----------|-------------|
| 1 | Label Analysis | Inspect current segmentation state |
| 2 | Remove Isolated Voxels | Delete small isolated components |
| 3 | Remove Low Intensity | Delete voxels with HU ≤ 0 |
| 4 | Remove High Intensity | Delete voxels with HU ≥ threshold |
| 5 | Smoothing | Smooth label boundaries |
| 6 | Boundary Expansion | Expand boundaries based on intensity |
| 7 | Boundary Trimming | Trim surface voxels outside HU range |
| 8 | Staircase Fill | Fill staircase-like boundary gaps |
| 9 | Protrusion Removal | Remove thin protrusions |
| 10 | Fill Internal Holes | Fill holes inside labels |
| 11 | Relabel Isolated Kidney → Tumor | Convert tumor-adjacent isolated kidney to tumor |
| 12 | Convex Hull Labeling | Create 3D convex hull from seeds and fill interior |
| 13 | Merge Segmentations | Replace kidney label from external file |
| 14 | Phase Comparison | Compare all phases simultaneously |
| r | Rollback | Revert to previous state |

---

## Detailed Function Descriptions

### Function 1: Label Analysis

Displays the current state of the segmentation. Does not modify data.

**Output:**
- Voxel count per label
- Intensity mean ± std (HU) per label (if CT image is available)
- Connected component count for kidney/tumor
- Size of top 5 components
- Surface ratio — higher values indicate more irregular boundaries

**No additional input required** — runs immediately.

---

### Function 2: Remove Isolated Voxels

Analyzes connected components of kidney or tumor labels, keeps only the top N largest, and removes the rest to background (0).

**Choices:**

| Input | Target | Components Kept |
|-------|--------|----------------|
| `1` | Kidney (label 1) | Top 2 (left/right kidney) |
| `2` | Tumor (label 2) | Top 1 |
| `3` | Both | Kidney top 2 + Tumor top 1 |

**No additional parameters** — runs after selection.

---

### Function 3: Remove Low Intensity

Removes kidney/tumor voxels with CT intensity ≤ 0 HU (air, fat, non-organ tissue).

**Requires:** CT image file must be available.

**No additional input** — runs immediately with fixed threshold of 0.

---

### Function 4: Remove High Intensity

Removes kidney/tumor voxels with CT intensity ≥ threshold (calcification, contrast-enhanced vessels).

**Requires:** CT image file must be available.

**Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Threshold | Voxels with HU ≥ this value are removed | 400 |

---

### Function 5: Smoothing

Combines morphological operations (closing, opening) with Gaussian smoothing to create smoother boundaries.

**Step 1 — Target Selection:**

| Input | Target | Behavior |
|-------|--------|----------|
| `1` | Tumor (label 2) | Smooths tumor boundary only. Keeps only the largest component. Constrained within kidney+tumor region |
| `2` | Cyst (label 3) | Smooths cyst boundary. Constrained within kidney+cyst region. Auto-fills internal holes |
| `3` | Whole organ surface | Merges kidney+tumor+cyst and smooths the outer surface. Preserves tumor/cyst labels inside |

**Step 2 — Parameter Input:**

| Parameter | Description | Tumor Default | Cyst Default | Organ Default |
|-----------|-------------|--------------|-------------|---------------|
| Gaussian sigma (mm) | Smoothing strength. Larger = smoother | 1.0 | 1.0 | 1.0 |
| Closing iterations | Gap-filling strength | 3 | 2 | 3 |
| Opening iterations | Protrusion removal strength | 2 | 1 | 2 |
| Components to keep | (Organ only) Keep top N components | - | - | 2 |

**Tuning Tips:**
- Rough boundary → increase sigma (1.5–2.0)
- Many gaps → increase closing iterations (4–5)
- Over-eroded / concave → decrease opening iterations (0–1)
- Excessive volume loss → decrease opening iterations

---

### Function 6: Boundary Expansion

Iteratively dilates the selected label boundary by 1 voxel per step, absorbing only voxels that meet the intensity condition.

**Step 1 — Target Selection:**

| Input | Target | Expandable Region |
|-------|--------|-------------------|
| `1` | Kidney (label 1) | Background (0) only |
| `2` | Tumor (label 2) | Background (0) + Kidney (1) |
| `3` | Cyst (label 3) | Background (0) + Kidney (1) |

After selection, the current label's intensity statistics (mean ± std HU) are displayed automatically.

**Step 2 — Intensity Condition:**

| Input | Condition | Use Case |
|-------|-----------|----------|
| `1` | Lower bound only (≥ threshold) | Absorb voxels above a minimum HU. Good for kidney expansion |
| `2` | Bidirectional range (mean ± tolerance) | Absorb voxels with similar HU. Good for cyst/tumor expansion |

**Step 3 — Parameter Input:**

Lower bound mode:

| Parameter | Description | Default |
|-----------|-------------|---------|
| Minimum intensity HU | Only absorb voxels ≥ this value | 120 |

Bidirectional range mode:

| Parameter | Description | Default |
|-----------|-------------|---------|
| HU range (mean ± X) | Only absorb voxels within this range of the mean | max(std×2, 15) |

Common:

| Parameter | Description | Default |
|-----------|-------------|---------|
| Expansion steps | Repeat dilation this many times (1 voxel per step) | 5 |

---

### Function 7: Boundary Trimming

Trims surface voxels that fall outside the HU range, working inward from the outer surface one layer at a time.

**Step 1 — Target Selection:**

| Input | Target | Trimmed Voxel Handling |
|-------|--------|----------------------|
| `1` | Kidney — whole organ surface (kidney+tumor+cyst) | All become background (0) |
| `2` | Tumor (label 2) | Adjacent to background → background (0), otherwise → kidney (1) |
| `3` | Cyst (label 3) | Adjacent to background → background (0), otherwise → kidney (1) |

After selection, the target's intensity statistics (mean ± std HU) are displayed automatically.

**Step 2 — Parameter Input:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| HU range (mean ± X) | Surface voxels outside this range are trimmed | max(std×2, 15) |
| Max iterations | Number of surface layers to trim | 1 |

---

### Function 8: Staircase Fill

Fills staircase-like gaps in the organ boundary (kidney+tumor) using 26-connectivity binary closing. Filled voxels become kidney (1). Tumor and cyst labels are protected.

**Requires:** HU filter is applied when CT image is available.

**Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Closing iterations | Fill strength. 1 fills single-voxel diagonal gaps | 1 |
| Minimum intensity HU | Only fill voxels ≥ this value | 120 |

---

### Function 9: Protrusion Removal

Removes thin protrusions from the kidney (label 1) boundary using morphological opening. Kidney voxels adjacent to tumor are protected.

**Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Opening iterations | Removal strength. Higher removes thicker protrusions | 1 |

---

### Function 10: Fill Internal Holes

Fills holes inside labels using `binary_fill_holes`.

**Choices:**

| Input | Target | Behavior |
|-------|--------|----------|
| `1` | Kidney (label 1) | Fill kidney holes as kidney |
| `2` | Tumor (label 2) | Fill tumor holes as tumor |
| `3` | Kidney+Tumor whole organ | Merge and fill. Filled voxels become kidney (1). Cyst protected |
| `4` | Cyst (label 3) | Fill cyst holes as cyst |

**No additional parameters** — runs after selection.

---

### Function 11: Relabel Isolated Kidney → Tumor

Finds kidney connected components beyond the top 2 (left/right kidney body). Among these isolated components, those adjacent to tumor are relabeled as tumor (2).

Used to correct cases where tumor was mistakenly labeled as kidney.

**No additional input** — runs immediately.

---

### Function 12: Convex Hull Labeling

Generates a 3D Convex Hull from seed voxels pre-painted in an image viewer, then fills the interior.

**Preparation:** Use an image viewer (e.g., ITK-SNAP) to paint seeds — label 2 for tumor, label 3 for cyst — across multiple planes (axial/sagittal/coronal).

**Step 1 — Target Selection:**

| Input | Target | Seed Label | Allowed Region | Protected Label |
|-------|--------|-----------|----------------|-----------------|
| `1` | Tumor (label 2) | label 2 | Kidney + Tumor | Cyst (3) protected |
| `2` | Cyst (label 3) | label 3 | Kidney + Cyst | Tumor (2) protected |

**No additional parameters** — Convex Hull computation and filling runs automatically after selection.

**Note:** If seeds exist on only one plane, the convex hull will be degenerate and fail. Always paint seeds across multiple planes.

---

### Function 13: Merge Segmentations

Imports kidney (label 1) from an external NIfTI file and replaces the current segmentation's kidney. Tumor (2) and cyst (3) from the current file are preserved.

Used to combine kidney/tumor segmentations from different models or methods.

**Input:**

| Parameter | Description |
|-----------|-------------|
| File path | Path to the NIfTI file to import kidney label from |

**Note:** The external file's shape must match the current file.

---

### Function 14: Phase Comparison

Loads all phases (A, D, P) simultaneously and performs cross-phase analysis. Does not modify data. Runs regardless of phase selection.

**Output:**

1. **Label Size Comparison** — Voxel counts for kidney/tumor/cyst/total organ per phase, displayed side by side. Maximum difference (%) between phases is shown.

2. **Pairwise Dice Coefficients** — Dice score (0–1) for each label across all phase pairs (A↔D, A↔P, D↔P). Closer to 1 means better agreement.

3. **Disagreement Slice Analysis** — For each phase pair and label, shows the range of slices with differences and the top 5 slices with the largest disagreement.

**No additional input** — runs immediately.

---

### Rollback (r)

Reverts to the previous state. Multiple rollback levels are supported.
After any function execution, you can immediately enter `r` to rollback.

---

## Dependencies

```bash
pip install nibabel numpy scipy
```
