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

- Enter `b` at any input step to cancel and return to the function selection menu.
- Rollback (`r`) is available at the function selection menu or immediately after execution.
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
| 14 | Resample & Merge | Resample external segmentation with different shape and merge |
| 15 | Phase Comparison | Compare all phases simultaneously |
| m | Region Restricted Run | Restrict to slice/bbox/direction then run a function |
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
| `1` | Tumor (label 2) | Smooths tumor boundary only. Keeps only the largest component. Cyst protected. Lost voxels assigned to background/kidney based on proximity |
| `2` | Cyst (label 3) | Smooths cyst boundary. Constrained within kidney+cyst region. Auto-fills internal holes. Lost voxels assigned to background/kidney based on proximity |
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

After selection, the target's intensity statistics (mean ± std, min, max HU) are displayed automatically.

**Step 2 — Intensity Condition:**

| Input | Condition | Description |
|-------|-----------|-------------|
| `1` | Range | Remove outside mean ± tolerance (default: max(std×2, 15)) |
| `2` | Lower bound | Remove HU < threshold (default: 0) |
| `3` | Upper bound | Remove HU > threshold (default: 400) |

**Step 3 — Iterations:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Max iterations | Number of surface layers to trim | 1 |

---

### Function 8: Staircase Fill

Fills staircase-like gaps in the organ boundary (kidney+tumor+cyst).
Filled voxels are assigned the nearest existing label (kidney/tumor/cyst).

**Requires:** HU filter is applied when CT image is available.

**Step 1 — Method Selection:**

| Input | Method | Description |
|-------|--------|-------------|
| `1` | Closing | 26-connectivity closing to fill concave gaps. 1 iteration fills single-voxel diagonal gaps |
| `2` | Convex Hull | Per-slice 2D Convex Hull to fill convex staircase patterns |

**Common Parameter:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Minimum intensity HU | Only fill voxels ≥ this value | 120 |

**Closing mode additional parameter:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| Closing iterations | Fill strength | 1 |

**Convex Hull mode additional parameter:**

| Parameter | Description |
|-----------|-------------|
| Slice axis | Select axis 0 / 1 / 2 |

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

Generates a Convex Hull from seed voxels pre-painted in an image viewer, then fills the interior.

**Preparation:** Use an image viewer (e.g., 3D Slicer, ITK-SNAP) to paint seeds — label 2 for tumor, label 3 for cyst.

**Step 1 — Target Selection:**

| Input | Target | Seed Label | Protected Label |
|-------|--------|-----------|-----------------|
| `1` | Tumor (label 2) | label 2 | Cyst (3) protected |
| `2` | Cyst (label 3) | label 3 | Tumor (2) protected |

All background/kidney voxels inside the hull are filled. Only the protected label is preserved.

**Step 1.5 — Component Selection (when 2+ components exist):**

If seeds are split into multiple disconnected components, a list is displayed:

```
Cyst component 3 found:
  1: size 1,234 voxels (slice 45~78)
  2: size 567 voxels (slice 90~110)
  3: size 89 voxels (slice 120~125)

Select:
  a: Process all individually
  1,2,3: Select by number (comma-separated for multiple)
```

- Enter `a`: Each component gets its own independent convex hull
- Enter numbers: Only selected components are processed (e.g., `1,3` for 1st and 3rd)

This step is skipped when there is only 1 component.

**Step 2 — Method Selection:**

| Input | Method | Best For |
|-------|--------|----------|
| `1` | 3D ConvexHull | Seeds painted across multiple axes (axial/sagittal/coronal). Builds a 3D convex hull |
| `2` | Slice-by-slice 2D ConvexHull + interpolation | Seeds painted along a single axis across multiple slices. Builds 2D hulls per slice and interpolates between them |

**2D method additional input:**

| Parameter | Description |
|-----------|-------------|
| Slice axis | Select the axis along which seeds were painted (axis 0 / 1 / 2) |

**Tips:**
- 3D method: Seeds must be well-distributed across multiple axes. If seeds lie on a single plane, the hull degenerates and produces 0 results
- 2D method: Requires at least 2 seed slices for interpolation to work. Even widely spaced slices are automatically filled between

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

### Function 14: Resample & Merge

Resamples labels from an external segmentation file with a different shape into the current segmentation's coordinate space. Useful when phases have different shapes (e.g., applying kidney labels from A phase to D/P phase).

Function 13 (Merge Segmentations) requires matching shapes, while this function uses the NIfTI affine matrix to perform spatial coordinate transformation with nearest neighbor resampling.

**Step 1 — Source file path input:**

Source and current file shape/voxel sizes are displayed.

**Step 2 — Label Selection:**

| Input | Target | Behavior |
|-------|--------|----------|
| `1` | Kidney (label 1) | Replace current kidney with source kidney |
| `2` | Tumor (label 2) | Replace current tumor with source tumor |
| `3` | Cyst (label 3) | Replace current cyst with source cyst |
| `4` | All (kidney+tumor+cyst) | Import all labels from source |

Labels not selected for import are preserved from the current file.

---

### Function 15: Phase Comparison

Loads all phases (A, D, P) simultaneously and performs cross-phase analysis. Does not modify data. Runs regardless of phase selection.

**Output:**

1. **Label Size Comparison** — Voxel counts for kidney/tumor/cyst/total organ per phase, displayed side by side. Maximum difference (%) between phases is shown.

2. **Pairwise Dice Coefficients** — Dice score (0–1) for each label across all phase pairs (A↔D, A↔P, D↔P). Closer to 1 means better agreement.

3. **Disagreement Slice Analysis** — For each phase pair and label, shows the range of slices with differences and the top 5 slices with the largest disagreement.

**No additional input** — runs immediately.

---

### Region Restricted Run (m)

Restricts function execution to a specific region. Areas outside the region are preserved unchanged.

**Region restriction methods (can be combined):**

**1. Slice Range**

| Parameter | Description |
|-----------|-------------|
| Axis | Select axis 0 / 1 / 2 |
| Start slice | Range start |
| End slice | Range end |

**2. 3D Bounding Box**

| Parameter | Description |
|-----------|-------------|
| X (axis 0) start/end | X range |
| Y (axis 1) start/end | Y range |
| Z (axis 2) start/end | Z range |

**3. Direction Restriction**

Anatomical directions (R/L/A/P/S/I) are automatically detected from the NIfTI affine matrix, matching the direction labels shown in CT viewers.

| Direction | Meaning |
|-----------|---------|
| R (Right) / L (Left) | Lateral |
| A (Anterior) / P (Posterior) | Front/Back |
| S (Superior) / I (Inferior) | Up/Down |

Select one of 6 directions and a cut slice to restrict the region to that half.

**After region setup**, select a function (2–13) to run within the restricted region.

---

### Rollback (r)

Reverts to the previous state. Multiple rollback levels are supported.
Enter `r` at the function selection menu or immediately after execution to rollback.

---

## Dependencies

```bash
pip install nibabel numpy scipy
```
