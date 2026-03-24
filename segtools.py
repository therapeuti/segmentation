"""
세그멘테이션 후처리 대화형 도구

사용법:
    python segtools.py S004

    케이스 폴더를 지정하면 phase 선택 → 기능 선택 → 실행 → 반복
"""

import os
import sys
import re
import glob
import shutil
import nibabel as nib
import numpy as np
from scipy import ndimage
from scipy.spatial import ConvexHull, Delaunay


# ──────────────────────────────────────────────
# 유틸리티
# ──────────────────────────────────────────────

def surface_ratio(binary_mask):
    total = int(np.sum(binary_mask))
    if total == 0:
        return 0.0
    eroded = ndimage.binary_erosion(binary_mask)
    surface = total - int(np.sum(eroded))
    return surface / total * 100


def load_case(case_dir):
    """케이스 폴더에서 phase별 파일 경로 매칭."""
    case_name = os.path.basename(case_dir)
    phases = {}
    for phase in ["A", "D", "P"]:
        # 대소문자 둘 다 시도
        for fmt in [f"{case_name}_Segmentation_{phase}.nii.gz",
                    f"{case_name}_segmentation_{phase}.nii.gz"]:
            seg_path = os.path.join(case_dir, fmt)
            if os.path.exists(seg_path):
                img_path = os.path.join(case_dir, f"{case_name}_image_{phase}.nii.gz")
                phases[phase] = {
                    "seg": seg_path,
                    "img": img_path if os.path.exists(img_path) else None,
                }
                break
    return phases


def backup_file(path):
    """백업 생성."""
    backup_dir = os.path.join(os.path.dirname(os.path.abspath(path)), "backup_original")
    os.makedirs(backup_dir, exist_ok=True)
    backup_path = os.path.join(backup_dir, os.path.basename(path))
    if not os.path.exists(backup_path):
        shutil.copy2(path, backup_path)
        print(f"  Backup: {backup_path}")


def save_result(path, data, img):
    """결과 저장."""
    new_header = img.header.copy()
    new_header.set_data_dtype(np.uint16)
    new_header["scl_slope"] = 1
    new_header["scl_inter"] = 0
    new_img = nib.Nifti1Image(data, img.affine, new_header)
    nib.save(new_img, path)
    print(f"  저장 완료 Saved: {path}")


# phase별 롤백 히스토리 {phase: [(data, description), ...]}
rollback_history = {}



def print_label_info(data, ct_data=None):
    """현재 라벨 상태 출력."""
    labels = np.unique(data)
    print(f"\n  현재 라벨 상태 Current label status:")
    for label in labels:
        count = int(np.sum(data == label))
        name = {0: "배경 BG", 1: "신장 Kidney", 2: "종양 Tumor", 3: "물혹 Cyst"}.get(int(label), f"label{label}")
        info = f"    {name}(label {int(label)}): {count:>12,} voxels"
        if ct_data is not None and label > 0:
            vals = ct_data[data == label]
            info += f"  | intensity: {np.mean(vals):.0f} ± {np.std(vals):.0f} HU"
        print(info)
    print()


# ──────────────────────────────────────────────
# 기능 1: 고립 복셀 제거
# ──────────────────────────────────────────────

def func_remove_isolated(data, **kwargs):
    """고립 복셀 제거 (신장/종양 선택)"""
    target = input_choice("  대상 선택 Select target", ["1: 신장 Kidney (label 1)", "2: 종양 Tumor (label 2)", "3: 둘 다 Both"])
    keep_n_map = {"1": (1, 2), "2": (2, 1), "3": None}

    if target == "3":
        # 신장 먼저
        data = _remove_isolated_label(data, label=1, keep_n=2)
        data = _remove_isolated_label(data, label=2, keep_n=1)
    elif target == "1":
        data = _remove_isolated_label(data, label=1, keep_n=2)
    elif target == "2":
        data = _remove_isolated_label(data, label=2, keep_n=1)

    return data


def _remove_isolated_label(data, label, keep_n):
    """특정 라벨의 상위 N개 component만 유지."""
    name = {1: "신장", 2: "종양"}.get(label, f"label{label}")
    mask = (data == label)
    total = int(np.sum(mask))
    if total == 0:
        print(f"  {name}: 라벨 없음")
        return data

    labeled, n_comp = ndimage.label(mask)
    if n_comp <= keep_n:
        print(f"  {name}: {n_comp}개 component — 제거 대상 없음")
        return data

    sizes = ndimage.sum(mask, labeled, range(1, n_comp + 1))
    top_indices = np.argsort(sizes)[::-1][:keep_n]
    top_labels = set(top_indices + 1)

    result = data.copy()
    removed = 0
    for i in range(1, n_comp + 1):
        if i not in top_labels:
            comp_mask = (labeled == i)
            removed += int(np.sum(comp_mask))
            result[comp_mask] = 0

    print(f"  {name}: {n_comp}개 중 상위 {keep_n}개 유지, {n_comp - keep_n}개 제거 ({removed:,} voxels)")
    return result


# ──────────────────────────────────────────────
# 기능 2: 저강도 제거 (intensity ≤ 0)
# ──────────────────────────────────────────────

def func_remove_low_intensity(data, ct_data=None, **kwargs):
    """intensity ≤ 0인 신장/종양 복셀 삭제."""
    if ct_data is None:
        print("  CT 이미지 없음 No CT image — 실행 불가 cannot run")
        return data

    low = (ct_data <= 0)
    result = data.copy()

    for label, name in [(1, "신장"), (2, "종양")]:
        mask = (data == label) & low
        removed = int(np.sum(mask))
        result[mask] = 0
        if removed > 0:
            print(f"  {name}: -{removed:,} voxels (intensity ≤ 0)")
        else:
            print(f"  {name}: 해당 없음")

    return result


# ──────────────────────────────────────────────
# 기능 3: 고강도 제거 (intensity ≥ 400)
# ──────────────────────────────────────────────

def func_remove_high_intensity(data, ct_data=None, **kwargs):
    """intensity ≥ 400인 신장/종양 복셀 삭제."""
    if ct_data is None:
        print("  CT 이미지 없음 No CT image — 실행 불가 cannot run")
        return data

    threshold = input_int("  Threshold (기본 default 400)", default=400)
    high = (ct_data >= threshold)
    result = data.copy()

    for label, name in [(1, "신장"), (2, "종양")]:
        mask = (data == label) & high
        removed = int(np.sum(mask))
        result[mask] = 0
        if removed > 0:
            print(f"  {name}: -{removed:,} voxels (intensity ≥ {threshold})")
        else:
            print(f"  {name}: 해당 없음")

    return result


# ──────────────────────────────────────────────
# 기능: Smoothing (종양 / 물혹 / 장기 전체 외곽)
# ──────────────────────────────────────────────

def func_smooth(data, zooms=None, **kwargs):
    """Smoothing 통합 — 대상 선택 후 실행."""
    target = input_choice("  Smoothing 대상 Select target", [
        "1: 종양 Tumor (label 2)",
        "2: 물혹 Cyst (label 3)",
        "3: 장기 전체 외곽 Whole organ surface",
    ])

    if target == "1":
        return _smooth_tumor(data, zooms)
    elif target == "2":
        return _smooth_cyst(data, zooms)
    else:
        return _smooth_organ(data, zooms)


def _smooth_tumor(data, zooms):
    """종양 라벨 smoothing."""
    kidney_mask = (data == 1)
    tumor_mask = (data == 2)
    organ_mask = kidney_mask | tumor_mask
    before = int(np.sum(tumor_mask))

    if before == 0:
        print("  종양 라벨 없음")
        return data

    sigma = input_float("  Gaussian sigma mm (기본 default 1.0)", default=1.0)
    close_iter = input_int("  Closing 반복 iterations (기본 default 3)", default=3)
    open_iter = input_int("  Opening 반복 iterations (기본 default 2)", default=2)

    mask = tumor_mask.astype(np.float64)

    # 가장 큰 component만 유지
    labeled, n_comp = ndimage.label(mask)
    if n_comp > 1:
        sizes = ndimage.sum(mask, labeled, range(1, n_comp + 1))
        largest = np.argmax(sizes) + 1
        mask = (labeled == largest).astype(np.float64)
        print(f"  고립 제거: {n_comp - 1}개 component 삭제")

    struct = ndimage.generate_binary_structure(3, 1)
    if close_iter > 0:
        mask = ndimage.binary_closing(mask, structure=struct, iterations=close_iter).astype(np.float64)
    if open_iter > 0:
        mask = ndimage.binary_opening(mask, structure=struct, iterations=open_iter).astype(np.float64)

    sigma_voxels = [sigma / float(z) for z in zooms]
    smoothed = ndimage.gaussian_filter(mask, sigma=sigma_voxels)
    mask_final = (smoothed >= 0.5)

    result = data.copy()
    # 기존 종양 중 smoothed 결과에서 빠진 복셀 → 원래 라벨 복원
    lost = tumor_mask & ~mask_final
    result[lost] = 0  # 일단 배경으로
    # lost 중 신장에 둘러싸인 복셀만 신장으로
    if int(np.sum(lost)) > 0:
        from scipy.ndimage import distance_transform_edt
        bg_dist = distance_transform_edt(data != 0)
        kid_dist = distance_transform_edt(data != 1)
        result[lost & (kid_dist < bg_dist)] = 1
    # smoothed 종양 적용 — 배경 확장 허용, 다른 라벨(물혹 등)은 보호
    new_tumor = mask_final & (result != 3)
    result[new_tumor] = 2

    after = int(np.sum(result == 2))
    before_s = surface_ratio(tumor_mask.astype(np.uint8))
    after_s = surface_ratio((result == 2).astype(np.uint8))
    print(f"  Voxels: {before:,} -> {after:,} ({(after-before)/max(before,1)*100:+.1f}%)")
    print(f"  Surface: {before_s:.1f}% -> {after_s:.1f}%")

    return result


def _smooth_cyst(data, zooms):
    """물혹(label 3) 경계를 스무스하게 정리."""
    cyst_mask = (data == 3)
    kidney_mask = (data == 1)
    before = int(np.sum(cyst_mask))

    if before == 0:
        print("  물혹 라벨 없음")
        return data

    sigma = input_float("  Gaussian sigma mm (기본 default 1.0)", default=1.0)
    close_iter = input_int("  Closing 반복 iterations (기본 default 2)", default=2)
    open_iter = input_int("  Opening 반복 iterations (기본 default 1)", default=1)

    # 물혹이 존재할 수 있는 영역: 기존 신장 + 물혹
    allowed = kidney_mask | cyst_mask

    mask = cyst_mask.astype(np.float64)
    struct = ndimage.generate_binary_structure(3, 1)

    if close_iter > 0:
        mask = ndimage.binary_closing(mask, structure=struct, iterations=close_iter).astype(np.float64)
    if open_iter > 0:
        mask = ndimage.binary_opening(mask, structure=struct, iterations=open_iter).astype(np.float64)

    if zooms is not None:
        sigma_voxels = [sigma / float(z) for z in zooms]
    else:
        sigma_voxels = sigma

    smoothed = ndimage.gaussian_filter(mask, sigma=sigma_voxels)
    mask_final = (smoothed >= 0.5) & allowed

    # 내부 구멍 채우기
    mask_final = ndimage.binary_fill_holes(mask_final) & allowed

    result = data.copy()
    # 기존 물혹 중 smoothed 결과에서 빠진 복셀 처리
    lost = cyst_mask & ~mask_final
    if int(np.sum(lost)) > 0:
        # 원본에서 물혹을 제외한 라벨 기준으로 nearest label 할당
        from scipy.ndimage import distance_transform_edt
        for lbl, lbl_name in [(0, "배경"), (1, "신장")]:
            lbl_mask = (data == lbl)
            if not np.any(lbl_mask):
                continue
            if lbl == 0:
                bg_dist = distance_transform_edt(~lbl_mask)
            else:
                kid_dist = distance_transform_edt(~lbl_mask)
        # 배경과 신장 중 가까운 쪽으로
        use_bg = bg_dist <= kid_dist
        result[lost & use_bg] = 0
        result[lost & ~use_bg] = 1
    result[mask_final] = 3

    after = int(np.sum(result == 3))
    before_s = surface_ratio(cyst_mask.astype(np.uint8))
    after_s = surface_ratio((result == 3).astype(np.uint8))
    print(f"  Voxels: {before:,} → {after:,} ({(after - before) / max(before, 1) * 100:+.1f}%)")
    print(f"  Surface: {before_s:.1f}% → {after_s:.1f}%")

    return result


def _smooth_organ(data, zooms):
    """신장+종양+물혹 장기 전체 외곽 smoothing."""
    kidney_mask = (data == 1)
    tumor_mask = (data == 2)
    before_kidney = int(np.sum(kidney_mask))

    if before_kidney == 0:
        print("  신장 라벨 없음")
        return data

    sigma = input_float("  Gaussian sigma mm (기본 default 1.0)", default=1.0)
    close_iter = input_int("  Closing 반복 iterations (기본 default 3)", default=3)
    open_iter = input_int("  Opening 반복 iterations (기본 default 2)", default=2)
    keep_n = input_int("  유지할 component 수 Components to keep (기본 default 2)", default=2)

    cyst_mask = (data == 3)
    organ_mask = (kidney_mask | tumor_mask | cyst_mask).astype(np.uint8)

    # 상위 keep_n개 component 유지
    labeled, n_comp = ndimage.label(organ_mask)
    if n_comp > keep_n:
        sizes = ndimage.sum(organ_mask, labeled, range(1, n_comp + 1))
        top_indices = np.argsort(sizes)[::-1][:keep_n]
        top_labels = set(top_indices + 1)
        new_mask = np.zeros_like(organ_mask)
        for lbl in top_labels:
            new_mask[labeled == lbl] = 1
        removed = int(np.sum(organ_mask)) - int(np.sum(new_mask))
        organ_mask = new_mask.astype(np.uint8)
        labeled, n_comp = ndimage.label(organ_mask)
        print(f"  고립 제거: {removed:,} voxels")

    # 각 component별 smoothing
    struct = ndimage.generate_binary_structure(3, 1)
    smoothed_organ = np.zeros_like(organ_mask, dtype=np.uint8)

    for i in range(1, n_comp + 1):
        comp = (labeled == i).astype(np.float64)
        if close_iter > 0:
            comp = ndimage.binary_closing(comp, structure=struct, iterations=close_iter).astype(np.float64)
        if open_iter > 0:
            comp = ndimage.binary_opening(comp, structure=struct, iterations=open_iter).astype(np.float64)
        sigma_voxels = [sigma / float(z) for z in zooms]
        comp = ndimage.gaussian_filter(comp, sigma=sigma_voxels)
        smoothed_organ = np.maximum(smoothed_organ, (comp >= 0.5).astype(np.uint8))

    result = data.copy()
    result[kidney_mask] = 0
    result[tumor_mask] = 0
    result[cyst_mask] = 0
    final_tumor = tumor_mask & (smoothed_organ == 1)
    final_cyst = cyst_mask & (smoothed_organ == 1)
    final_kidney = (smoothed_organ == 1) & ~final_tumor & ~final_cyst
    result[final_kidney] = 1
    result[final_tumor] = 2
    result[final_cyst] = 3

    after_kidney = int(np.sum(result == 1))
    before_s = surface_ratio(organ_mask)
    after_s = surface_ratio((smoothed_organ == 1).astype(np.uint8))
    print(f"  신장 Voxels: {before_kidney:,} -> {after_kidney:,} ({(after_kidney-before_kidney)/max(before_kidney,1)*100:+.1f}%)")
    print(f"  장기 외곽 Surface: {before_s:.1f}% -> {after_s:.1f}%")

    return result


# ──────────────────────────────────────────────
# 기능: Intensity 기반 경계 확장 (신장/종양/물혹)
# ──────────────────────────────────────────────

def func_expand(data, ct_data=None, **kwargs):
    """Intensity 기반 경계 확장 — 대상 및 조건 선택."""
    if ct_data is None:
        print("  CT 이미지 없음 No CT image — 실행 불가 cannot run")
        return data

    target = input_choice("  확장 대상 Expansion target", [
        "1: 신장 Kidney (label 1)",
        "2: 종양 Tumor (label 2)",
        "3: 물혹 Cyst (label 3)",
    ])

    label = int(target)
    name = {1: "신장", 2: "종양", 3: "물혹"}[label]
    mask = (data == label)
    before = int(np.sum(mask))

    if before == 0:
        print(f"  {name} 라벨 없음")
        return data

    # 확장 가능 영역: 신장은 배경만, 종양/물혹은 배경+신장
    if label == 1:
        expandable = (data == 0)
        print("  확장 영역 Expandable: 배경 Background")
    else:
        expandable = (data == 0) | (data == 1)
        print("  확장 영역 Expandable: 배경 Background + 신장 Kidney")

    # intensity 통계 표시
    vals = ct_data[mask]
    val_mean = float(np.mean(vals))
    val_std = float(np.std(vals))
    print(f"  {name} intensity: {val_mean:.1f} ± {val_std:.1f} HU")

    # intensity 조건 선택
    mode = input_choice("  Intensity 조건 Condition", [
        "1: 하한만 Lower bound only (≥ threshold)",
        "2: 양방향 범위 Bidirectional range (mean ± tolerance)",
    ])

    if mode == "1":
        threshold = input_float("  최소 intensity HU Min intensity (기본 default 120)", default=120.0)
        hu_filter = (ct_data >= threshold)
        print(f"  조건: HU ≥ {threshold:.0f}")
    else:
        default_tol = round(max(val_std * 2, 15))
        tolerance = input_float(
            f"  허용 HU 범위 HU range (mean ± X, 기본 default {default_tol})",
            default=default_tol)
        hu_lo = val_mean - tolerance
        hu_hi = val_mean + tolerance
        hu_filter = (ct_data >= hu_lo) & (ct_data <= hu_hi)
        print(f"  조건: HU {hu_lo:.0f} ~ {hu_hi:.0f}")

    steps = input_int("  확장 횟수 Expansion steps (기본 default 5)", default=5)

    expand_mask = mask.copy()
    struct = ndimage.generate_binary_structure(3, 1)
    total_added = 0

    for step in range(steps):
        dilated = ndimage.binary_dilation(expand_mask, structure=struct)
        candidates = dilated & ~expand_mask & expandable
        accepted = candidates & hu_filter
        added = int(np.sum(accepted))
        total_added += added

        if added == 0:
            print(f"    Step {step + 1}: 확장 가능 복셀 없음 — 종료")
            break

        expand_mask[accepted] = True
        print(f"    Step {step + 1}: +{added:,} voxels")

    result = data.copy()
    new_voxels = expand_mask & ~mask
    result[new_voxels] = label
    print(f"  {name} 확장 완료: {before:,} → {before + total_added:,} (+{total_added:,} voxels)")

    return result


# ──────────────────────────────────────────────
# 기능 7: 라벨 상태 분석
# ──────────────────────────────────────────────

def func_analyze(data, ct_data=None, **kwargs):
    """현재 라벨 상태 상세 분석."""
    print_label_info(data, ct_data)

    for label, name in [(1, "신장"), (2, "종양")]:
        mask = (data == label)
        total = int(np.sum(mask))
        if total == 0:
            continue
        labeled, n_comp = ndimage.label(mask)
        sizes = ndimage.sum(mask, labeled, range(1, n_comp + 1))
        sorted_idx = np.argsort(sizes)[::-1]
        sr = surface_ratio(mask.astype(np.uint8))
        print(f"  {name}: {n_comp}개 component, surface={sr:.1f}%")
        for rank, idx in enumerate(sorted_idx[:5]):
            print(f"    #{rank+1}: {int(sizes[idx]):>10,} voxels")
        if n_comp > 5:
            small_total = sum(int(sizes[idx]) for idx in sorted_idx[5:])
            print(f"    ... 외 {n_comp-5}개 ({small_total:,} voxels)")
        print()

    return data  # 변경 없음


# ──────────────────────────────────────────────
# 입력 헬퍼
# ──────────────────────────────────────────────

class CancelOperation(Exception):
    """현재 기능 취소."""


def _check_special(val):
    """q/r/b 및 한글 대응(ㅂ/ㄱ/ㅠ) 입력 시 예외 발생."""
    if val.lower() in ("q", "b", "r", "ㅂ", "ㅠ", "ㄱ"):
        raise CancelOperation()


def input_choice(prompt, options):
    """선택지 입력."""
    for opt in options:
        print(f"    {opt}")
    print(f"    b: 돌아가기 Back")
    while True:
        val = input(f"  {prompt}: ").strip()
        _check_special(val)
        valid = [opt.split(":")[0].strip() for opt in options]
        if val in valid:
            return val
        print(f"    → {', '.join(valid)} 중 하나를 입력하세요")


def input_int(prompt, default):
    val = input(f"  {prompt} (b:돌아가기 Back): ").strip()
    _check_special(val)
    if val == "":
        return default
    try:
        return int(val)
    except ValueError:
        return default


def input_float(prompt, default):
    val = input(f"  {prompt} (b:돌아가기 Back): ").strip()
    _check_special(val)
    if val == "":
        return default
    try:
        return float(val)
    except ValueError:
        return default


# ──────────────────────────────────────────────
# 기능 8: 내부 구멍 채우기
# ──────────────────────────────────────────────

def func_fill_holes(data, **kwargs):
    """라벨 내부 구멍을 채우기 (외곽 형태 유지)."""
    target = input_choice("  대상 선택 Select target", ["1: 신장 Kidney (label 1)", "2: 종양 Tumor (label 2)", "3: 신장+종양 장기 전체 Whole organ", "4: 물혹 Cyst (label 3)"])

    if target == "3":
        # 신장+종양 합쳐서 fill → 채워진 부분은 신장으로
        kidney_mask = (data == 1)
        tumor_mask = (data == 2)
        organ_mask = (kidney_mask | tumor_mask)
        before = int(np.sum(organ_mask))

        filled = ndimage.binary_fill_holes(organ_mask).astype(np.uint8)
        new_voxels = (filled == 1) & ~organ_mask & (data != 3)  # 물혹 보호

        result = data.copy()
        result[new_voxels] = 1  # 채워진 부분은 신장으로
        added = int(np.sum(new_voxels))
        print(f"  장기 내부 구멍 채움: +{added:,} voxels (신장으로 라벨링)")
    else:
        label = int(target)
        name = {1: "신장", 2: "종양", 4: "물혹"}[label]
        if label == 4:
            label = 3
        mask = (data == label)
        before = int(np.sum(mask))

        filled = ndimage.binary_fill_holes(mask).astype(np.uint8)
        new_voxels = (filled == 1) & ~mask

        result = data.copy()
        result[new_voxels] = label
        added = int(np.sum(new_voxels))
        print(f"  {name} 내부 구멍 채움: +{added:,} voxels")

    return result


# ──────────────────────────────────────────────
# 기능 9: 고립 신장 → 종양 재라벨링
# ──────────────────────────────────────────────

def func_relabel_isolated_kidney(data, **kwargs):
    """종양 인접 고립 신장 component를 종양으로 재라벨링."""
    kidney_mask = (data == 1)
    tumor_mask = (data == 2)

    total_kidney = int(np.sum(kidney_mask))
    if total_kidney == 0:
        print("  신장 라벨 없음")
        return data

    labeled, n_comp = ndimage.label(kidney_mask)
    if n_comp <= 2:
        print(f"  신장 component {n_comp}개 — 고립 없음")
        return data

    sizes = ndimage.sum(kidney_mask, labeled, range(1, n_comp + 1))
    sorted_indices = np.argsort(sizes)[::-1]

    # 상위 2개 = 좌우 신장 본체
    main_labels = set()
    for idx in sorted_indices[:2]:
        main_labels.add(idx + 1)

    # 종양 인접 마스크
    struct = ndimage.generate_binary_structure(3, 1)
    tumor_adj = ndimage.binary_dilation(tumor_mask, structure=struct).astype(bool)

    result = data.copy()
    relabeled_count = 0
    relabeled_voxels = 0

    for i in range(1, n_comp + 1):
        if i in main_labels:
            continue
        comp_mask = (labeled == i)
        comp_size = int(np.sum(comp_mask))

        if np.any(comp_mask & tumor_adj):
            result[comp_mask] = 2
            relabeled_count += 1
            relabeled_voxels += comp_size

    if relabeled_voxels > 0:
        print(f"  고립 신장 {relabeled_count}개 → 종양으로 재라벨링 (+{relabeled_voxels:,} voxels)")
    else:
        print(f"  종양 인접 고립 신장 없음")

    return result


# ──────────────────────────────────────────────
# 기능: 볼록 껍질 기반 라벨링 (종양/물혹)
# ──────────────────────────────────────────────

def func_label_convex(data, ct_data=None, zooms=None, **kwargs):
    """3D 또는 슬라이스별 2D 볼록 껍질로 라벨링."""
    target = input_choice("  라벨링 대상 Labeling target", [
        "1: 종양 Tumor (label 2)",
        "2: 물혹 Cyst (label 3)",
    ])

    if target == "1":
        label = 2
        name = "종양"
        seed_mask = (data == 2)
        protect_mask = (data == 3)  # 물혹 보호
    else:
        label = 3
        name = "물혹"
        seed_mask = (data == 3)
        protect_mask = (data == 2)  # 종양 보호

    n_seed = int(np.sum(seed_mask))
    if n_seed == 0:
        print(f"  {name} 시드(label {label})가 없습니다.")
        print(f"  → 여러 평면(axial/sagittal/coronal)에서 {name} 영역을 label {label}로 칠해주세요.")
        return data

    # ── Component 분리 및 선택 ──
    labeled_arr, n_comp = ndimage.label(seed_mask)
    if n_comp >= 2:
        comp_info = []
        for ci in range(1, n_comp + 1):
            comp_mask = (labeled_arr == ci)
            comp_size = int(np.sum(comp_mask))
            comp_coords = np.argwhere(comp_mask)
            slices_range = f"slice {comp_coords[:, 0].min()}~{comp_coords[:, 0].max()}"
            comp_info.append((ci, comp_size, slices_range, comp_mask))

        comp_info.sort(key=lambda x: -x[1])  # 크기 순 정렬

        print(f"\n  {name} component {n_comp}개 발견:")
        for idx, (ci, sz, sr, _) in enumerate(comp_info):
            print(f"    {idx + 1}: 크기 size {sz:,} voxels ({sr})")

        options = [f"a: 전체 각각 자동 처리 Process all individually"]
        options += [f"{idx + 1}: component {idx + 1} ({comp_info[idx][1]:,} voxels)" for idx in range(len(comp_info))]
        sel = input_choice("  처리할 component Select component(s) (쉼표로 복수 선택 comma-separated)", options)

        if sel.lower() == 'a':
            selected_masks = [(f"component {idx+1}", ci[3]) for idx, ci in enumerate(comp_info)]
        else:
            indices = [int(s.strip()) - 1 for s in sel.split(",") if s.strip().isdigit()]
            indices = [i for i in indices if 0 <= i < len(comp_info)]
            if not indices:
                print("  유효한 선택 없음 No valid selection")
                return data
            selected_masks = [(f"component {i+1}", comp_info[i][3]) for i in indices]
    else:
        selected_masks = [("전체 all", seed_mask)]

    # ── 방식 선택 (한 번만) ──
    method = input_choice("  방식 Method", [
        "1: 3D ConvexHull (시드가 여러 축에 분포된 경우 When seeds span multiple axes)",
        "2: 슬라이스별 2D ConvexHull + 보간 Slice-by-slice 2D + interpolation",
    ])

    # 2D 방식이면 축도 한 번만 선택
    slice_axis = None
    if method == "2":
        axis_sel = input_choice("  슬라이스 축 Slice axis", [
            f"0: axis 0 (shape={data.shape[0]})",
            f"1: axis 1 (shape={data.shape[1]})",
            f"2: axis 2 (shape={data.shape[2]})",
        ])
        slice_axis = int(axis_sel)

    # ── 각 component별 처리 ──
    result = data.copy()
    total_added = 0

    for comp_name, comp_seed in selected_masks:
        comp_n = int(np.sum(comp_seed))
        print(f"\n  [{comp_name}] 시드 seeds: {comp_n:,} voxels")

        if method == "1":
            fill_mask = _label_convex_3d(data, comp_seed, comp_n)
        else:
            fill_mask = _label_convex_2d(data, comp_seed, comp_n, axis=slice_axis)

        if fill_mask is None:
            print(f"  [{comp_name}] 건너뜀 skipped")
            continue

        # 시드 포함 + 보호 라벨 침범 방지
        fill_mask = fill_mask | comp_seed
        fill_mask = (fill_mask & ~protect_mask) | comp_seed

        added = int(np.sum(fill_mask)) - comp_n
        total_added += added
        print(f"  [{comp_name}] {comp_n:,} → {comp_n + added:,} voxels (+{added:,})")

        result[fill_mask] = label

    print(f"\n  총 결과 Total: {n_seed:,} → {n_seed + total_added:,} voxels (+{total_added:,})")
    return result


def _label_convex_3d(data, seed_mask, n_seed):
    """3D ConvexHull로 내부 전체 채움."""
    coords = np.argwhere(seed_mask)
    print(f"  꼭짓점 Vertices: {len(coords):,} voxels")

    if len(coords) < 4:
        print("  시드 4개 미만 — 3D 볼록 껍질 불가")
        return None

    try:
        hull = ConvexHull(coords)
        delaunay = Delaunay(coords[hull.vertices])
    except Exception as e:
        print(f"  볼록 껍질 계산 실패 ConvexHull failed: {e}")
        print("  → 시드를 여러 축에 걸쳐 칠해주세요 Paint seeds across multiple axes")
        return None

    margin = 1
    mins = np.maximum(coords.min(axis=0) - margin, 0)
    maxs = np.minimum(coords.max(axis=0) + margin + 1, data.shape)

    grid = np.mgrid[mins[0]:maxs[0], mins[1]:maxs[1], mins[2]:maxs[2]]
    test_points = grid.reshape(3, -1).T

    print(f"  내부 판정 중 Testing interior... ({len(test_points):,} voxels)")
    inside = delaunay.find_simplex(test_points) >= 0

    fill_mask = np.zeros(data.shape, dtype=bool)
    fill_mask[test_points[inside, 0], test_points[inside, 1], test_points[inside, 2]] = True
    return fill_mask


def _label_convex_2d(data, seed_mask, n_seed, axis=None):
    """슬라이스별 2D ConvexHull + 슬라이스 간 보간."""
    if axis is None:
        axis = input_choice("  슬라이스 축 Slice axis", [
            f"0: axis 0 (shape={data.shape[0]})",
            f"1: axis 1 (shape={data.shape[1]})",
            f"2: axis 2 (shape={data.shape[2]})",
        ])
        axis = int(axis)

    # 시드가 있는 슬라이스 찾기
    seed_slices = []
    for i in range(data.shape[axis]):
        sl = [slice(None)] * 3
        sl[axis] = i
        if np.any(seed_mask[tuple(sl)]):
            seed_slices.append(i)

    if len(seed_slices) == 0:
        print("  시드 슬라이스 없음 No seed slices")
        return None

    print(f"  시드 슬라이스 Seed slices: {len(seed_slices)}개 (범위 range: {seed_slices[0]}~{seed_slices[-1]})")

    # 슬라이스별 2D Convex Hull
    hull_masks = {}
    other_axes = [s for i, s in enumerate(data.shape) if i != axis]

    for si in seed_slices:
        sl = [slice(None)] * 3
        sl[axis] = si
        slice_seed = seed_mask[tuple(sl)]
        coords_2d = np.argwhere(slice_seed)
        if len(coords_2d) < 3:
            hull_masks[si] = slice_seed.copy()
            continue
        try:
            hull = ConvexHull(coords_2d)
            delaunay = Delaunay(coords_2d[hull.vertices])
            grid = np.mgrid[0:other_axes[0], 0:other_axes[1]]
            test_pts = grid.reshape(2, -1).T
            inside = delaunay.find_simplex(test_pts) >= 0
            hull_mask = inside.reshape(other_axes[0], other_axes[1])
            hull_masks[si] = hull_mask | slice_seed
        except Exception:
            hull_masks[si] = slice_seed.copy()

    # 슬라이스 간 보간
    fill_mask = np.zeros(data.shape, dtype=bool)

    for si, hmask in hull_masks.items():
        sl = [slice(None)] * 3
        sl[axis] = si
        fill_mask[tuple(sl)] = hmask

    for idx in range(len(seed_slices) - 1):
        s_start = seed_slices[idx]
        s_end = seed_slices[idx + 1]
        if s_end - s_start <= 1:
            continue
        mask_start = hull_masks[s_start].astype(float)
        mask_end = hull_masks[s_end].astype(float)
        for si in range(s_start + 1, s_end):
            t = (si - s_start) / (s_end - s_start)
            interp = mask_start * (1 - t) + mask_end * t
            sl = [slice(None)] * 3
            sl[axis] = si
            fill_mask[tuple(sl)] = interp >= 0.5

    return fill_mask


# ──────────────────────────────────────────────
# 기능 11: 경계 계단 메꿈
# ──────────────────────────────────────────────

def func_fill_staircase(data, ct_data=None, **kwargs):
    """경계 계단 메꿈 — Closing 또는 슬라이스별 Convex Hull."""
    kidney_mask = (data == 1)
    tumor_mask = (data == 2)
    cyst_mask = (data == 3)
    full_organ = kidney_mask | tumor_mask | cyst_mask

    before = int(np.sum(full_organ))
    if before == 0:
        print("  장기 라벨 없음 No organ labels")
        return data

    mode = input_choice("  메꿈 방식 Fill method", [
        "1: Closing (오목한 틈새 채움 Fill concave gaps)",
        "2: Convex Hull (볼록하게 채움 Fill convex per slice)",
    ])

    threshold = input_float("  최소 intensity HU Min intensity (기본 default 120)", default=120.0)

    if mode == "1":
        return _fill_staircase_closing(data, ct_data, full_organ, threshold)
    else:
        return _fill_staircase_convex(data, ct_data, full_organ, threshold)


def _fill_staircase_closing(data, ct_data, full_organ, threshold):
    """Closing 기반 오목한 틈새 채움."""
    iterations = input_int("  Closing 반복 iterations (기본 default 1)", default=1)

    struct = ndimage.generate_binary_structure(3, 3)
    closed = ndimage.binary_closing(full_organ, structure=struct, iterations=iterations)

    new_voxels = closed & ~full_organ
    if ct_data is not None:
        new_voxels = new_voxels & (ct_data >= threshold)

    result = data.copy()
    if int(np.sum(new_voxels)) > 0:
        result = _assign_nearest_label(data, result, new_voxels)

    added = int(np.sum(new_voxels))
    added_kidney = int(np.sum(new_voxels & (result == 1)))
    added_tumor = int(np.sum(new_voxels & (result == 2)))
    added_cyst = int(np.sum(new_voxels & (result == 3)))
    print(f"  Closing 메꿈: +{added:,} voxels (신장 {added_kidney:,}, 종양 {added_tumor:,}, 물혹 {added_cyst:,}, HU ≥ {threshold:.0f})")
    return result


def _fill_staircase_convex(data, ct_data, full_organ, threshold):
    """슬라이스별 2D Convex Hull로 바깥쪽 계단을 볼록하게 채움."""
    axis = input_choice("  슬라이스 축 Slice axis", [
        f"0: axis 0 (shape={data.shape[0]})",
        f"1: axis 1 (shape={data.shape[1]})",
        f"2: axis 2 (shape={data.shape[2]})",
    ])
    axis = int(axis)

    convex_organ = np.zeros_like(full_organ, dtype=bool)

    n_slices = data.shape[axis]
    filled_slices = 0

    for i in range(n_slices):
        slicing = [slice(None)] * 3
        slicing[axis] = i
        sl = tuple(slicing)

        slice_mask = full_organ[sl]
        if not np.any(slice_mask):
            continue

        coords = np.argwhere(slice_mask)
        if len(coords) < 3:
            convex_organ[sl] = slice_mask
            continue

        try:
            hull = ConvexHull(coords)
            hull_path = Delaunay(coords[hull.vertices])
            # 슬라이스 전체 좌표에서 hull 내부 판정
            all_coords = np.argwhere(np.ones(slice_mask.shape, dtype=bool))
            inside = hull_path.find_simplex(all_coords) >= 0
            hull_mask = np.zeros(slice_mask.shape, dtype=bool)
            hull_mask[all_coords[inside, 0], all_coords[inside, 1]] = True
            convex_organ[sl] = hull_mask
            if np.any(hull_mask & ~slice_mask):
                filled_slices += 1
        except Exception:
            convex_organ[sl] = slice_mask

    new_voxels = convex_organ & ~full_organ
    if ct_data is not None:
        new_voxels = new_voxels & (ct_data >= threshold)

    result = data.copy()
    if int(np.sum(new_voxels)) > 0:
        result = _assign_nearest_label(data, result, new_voxels)

    added = int(np.sum(new_voxels))
    added_kidney = int(np.sum(new_voxels & (result == 1)))
    added_tumor = int(np.sum(new_voxels & (result == 2)))
    added_cyst = int(np.sum(new_voxels & (result == 3)))
    print(f"  Convex Hull 메꿈: +{added:,} voxels (신장 {added_kidney:,}, 종양 {added_tumor:,}, 물혹 {added_cyst:,})")
    print(f"  변경된 슬라이스: {filled_slices}개, HU ≥ {threshold:.0f}")
    return result


def _assign_nearest_label(data, result, new_voxels):
    """새 복셀에 가장 가까운 기존 라벨 할당."""
    from scipy.ndimage import distance_transform_edt
    min_dist = np.full(data.shape, np.inf)
    nearest_label = np.zeros(data.shape, dtype=np.uint16)
    for lbl in [1, 2, 3]:
        lbl_mask = (data == lbl)
        if not np.any(lbl_mask):
            continue
        dist = distance_transform_edt(~lbl_mask)
        closer = dist < min_dist
        min_dist[closer] = dist[closer]
        nearest_label[closer] = lbl
    result[new_voxels] = nearest_label[new_voxels]
    return result


# ──────────────────────────────────────────────
# 기능 12: 신장 돌출부 제거
# ──────────────────────────────────────────────

def func_remove_protrusion(data, **kwargs):
    """신장 경계에서 얇게 돌출된 부분을 제거."""
    kidney_mask = (data == 1)
    tumor_mask = (data == 2)

    before = int(np.sum(kidney_mask))
    if before == 0:
        print("  신장 라벨 없음")
        return data

    iterations = input_int("  Opening 반복 iterations (기본 default 1)", default=1)

    # Opening: erosion → dilation — 얇은 돌출부가 제거됨
    struct = ndimage.generate_binary_structure(3, 1)
    opened = ndimage.binary_opening(kidney_mask, structure=struct, iterations=iterations)

    # 종양에 인접한 신장은 보호 (종양 주변 1복셀)
    tumor_adj = ndimage.binary_dilation(tumor_mask, structure=struct)
    protected = kidney_mask & tumor_adj
    opened = opened | protected

    removed_mask = kidney_mask & ~opened
    removed = int(np.sum(removed_mask))

    result = data.copy()
    result[removed_mask] = 0

    print(f"  돌출부 제거: -{removed:,} voxels (신장 {before:,} → {before - removed:,})")
    return result


# ──────────────────────────────────────────────
# 기능: 경계 HU 트리밍 (신장(장기 외곽)/종양/물혹)
# ──────────────────────────────────────────────

def func_trim_boundary(data, ct_data=None, **kwargs):
    """경계에서 HU 범위 밖 복셀을 바깥부터 반복 깎아냄."""
    if ct_data is None:
        print("  CT 이미지 없음 No CT image — 실행 불가 cannot run")
        return data

    target = input_choice("  트리밍 대상 Trimming target", [
        "1: 장기 전체 외곽 Whole organ surface",
        "2: 종양 Tumor (label 2)",
        "3: 물혹 Cyst (label 3)",
    ])

    if target == "1":
        return _trim_organ(data, ct_data)
    elif target == "2":
        return _trim_single(data, ct_data, label=2, name="종양")
    else:
        return _trim_single(data, ct_data, label=3, name="물혹")


def _determine_removed_label(data, removed_mask):
    """깎인 복셀 각각의 주변을 확인하여 라벨 결정.

    배경(0)과 닿아있으면 배경으로, 아니면 신장(1)으로 변경.
    """
    struct = ndimage.generate_binary_structure(3, 1)
    bg_mask = (data == 0)
    bg_adj = ndimage.binary_dilation(bg_mask, structure=struct)
    # 배경과 닿아있는 깎인 복셀 → 배경(0), 나머지 → 신장(1)
    to_bg = removed_mask & bg_adj
    to_kidney = removed_mask & ~bg_adj
    return to_bg, to_kidney


def _trim_single(data, ct_data, label, name):
    """종양 또는 물혹 단일 라벨 트리밍."""
    mask = (data == label)
    before = int(np.sum(mask))
    if before == 0:
        print(f"  {name} 라벨 없음")
        return data

    vals = ct_data[mask]
    val_mean = float(np.mean(vals))
    val_std = float(np.std(vals))
    print(f"  {name}: {before:,} voxels, intensity {val_mean:.1f} ± {val_std:.1f} HU")

    default_tol = round(max(val_std * 2, 15))
    tolerance = input_float(
        f"  허용 HU 범위 HU range (mean ± X, 기본 default {default_tol})",
        default=default_tol)
    max_iter = input_int("  최대 반복 Max iterations (기본 default 1)", default=1)

    hu_lo = val_mean - tolerance
    hu_hi = val_mean + tolerance
    print(f"  HU 범위: {hu_lo:.0f} ~ {hu_hi:.0f}")

    hu_bad = (ct_data < hu_lo) | (ct_data > hu_hi)
    struct = ndimage.generate_binary_structure(3, 1)
    trimmed = mask.copy()

    for i in range(max_iter):
        eroded = ndimage.binary_erosion(trimmed, structure=struct)
        surface = trimmed & ~eroded
        to_remove = surface & hu_bad
        removed = int(np.sum(to_remove))
        if removed == 0:
            print(f"  {i + 1}회 반복 후 완료")
            break
        trimmed = trimmed & ~to_remove

    removed_mask = mask & ~trimmed
    after = int(np.sum(trimmed))
    print(f"  트리밍: {before:,} → {after:,} (−{before - after:,} voxels)")

    result = data.copy()
    to_bg, to_kidney = _determine_removed_label(data, removed_mask)
    result[to_bg] = 0
    result[to_kidney] = 1
    print(f"    → 배경으로: {int(np.sum(to_bg)):,}, 신장으로: {int(np.sum(to_kidney)):,}")

    return result


def _trim_organ(data, ct_data):
    """신장+종양+물혹 장기 전체 외곽 트리밍."""
    kidney_mask = (data == 1)
    tumor_mask = (data == 2)
    cyst_mask = (data == 3)
    organ_mask = kidney_mask | tumor_mask | cyst_mask
    before = int(np.sum(organ_mask))

    if before == 0:
        print("  장기 라벨 없음")
        return data

    vals = ct_data[organ_mask]
    val_mean = float(np.mean(vals))
    val_std = float(np.std(vals))
    print(f"  장기 전체: {before:,} voxels, intensity {val_mean:.1f} ± {val_std:.1f} HU")

    default_tol = round(max(val_std * 2, 15))
    tolerance = input_float(
        f"  허용 HU 범위 HU range (mean ± X, 기본 default {default_tol})",
        default=default_tol)
    max_iter = input_int("  최대 반복 Max iterations (기본 default 1)", default=1)

    hu_lo = val_mean - tolerance
    hu_hi = val_mean + tolerance
    print(f"  HU 범위: {hu_lo:.0f} ~ {hu_hi:.0f}")

    hu_bad = (ct_data < hu_lo) | (ct_data > hu_hi)
    struct = ndimage.generate_binary_structure(3, 1)
    trimmed = organ_mask.copy()

    for i in range(max_iter):
        eroded = ndimage.binary_erosion(trimmed, structure=struct)
        surface = trimmed & ~eroded
        # 장기 외곽만 깎음: 종양/물혹 내부는 보호
        # 표면 중 내부 라벨(종양/물혹) 경계는 제외하고 외곽만 대상
        to_remove = surface & hu_bad
        removed = int(np.sum(to_remove))
        if removed == 0:
            print(f"  {i + 1}회 반복 후 완료")
            break
        trimmed = trimmed & ~to_remove

    removed_mask = organ_mask & ~trimmed
    after = int(np.sum(trimmed))
    print(f"  트리밍: {before:,} → {after:,} (−{before - after:,} voxels)")

    # 깎인 복셀은 모두 배경으로 (장기 외곽이므로)
    result = data.copy()
    result[removed_mask] = 0
    print(f"    → 배경으로: {int(np.sum(removed_mask)):,}")

    return result


# ──────────────────────────────────────────────
# 기능: 세그멘테이션 합치기
# ──────────────────────────────────────────────

def func_merge_segmentations(data, **kwargs):
    """외부 파일에서 신장 라벨을 가져와 현재 세그멘테이션의 신장을 교체."""
    print("  신장(1) 라벨을 가져올 파일 경로 File path for kidney label:")
    path_a = input("  경로 Path: ").strip().strip('"')
    if not os.path.exists(path_a):
        print(f"  파일을 찾을 수 없음: {path_a}")
        raise CancelOperation()

    data_a = np.round(np.asanyarray(nib.load(path_a).dataobj)).astype(np.uint16)

    if data_a.shape != data.shape:
        print(f"  shape 불일치: 외부={data_a.shape}, 현재={data.shape}")
        raise CancelOperation()

    kidney_a = (data_a == 1)
    tumor_b = (data == 2)
    cyst_b = (data == 3)

    print(f"  외부 파일 — 신장: {int(np.sum(kidney_a)):,} voxels")
    print(f"  현재 파일 — 종양: {int(np.sum(tumor_b)):,}, 물혹: {int(np.sum(cyst_b)):,} voxels")

    # 현재 파일 기반으로 신장만 외부 파일 것으로 교체
    result = data.copy()
    result[result == 1] = 0       # 기존 신장 제거
    result[kidney_a & (result == 0)] = 1  # 외부 신장 삽입

    total = int(np.sum(result > 0))
    print(f"  합친 결과: 신장={int(np.sum(result == 1)):,}, "
          f"종양={int(np.sum(result == 2)):,}, "
          f"물혹={int(np.sum(result == 3)):,}, 전체={total:,} voxels")
    return result


# ──────────────────────────────────────────────
# 기능: Phase 비교 분석
# ──────────────────────────────────────────────

def func_compare_phases(phases):
    """모든 phase의 세그멘테이션을 로드하여 비교 분석."""
    phase_keys = sorted(phases.keys())
    if len(phase_keys) < 2:
        print("  비교할 phase가 2개 이상 필요합니다 Need at least 2 phases to compare.")
        return

    # ── 1. 모든 phase 로드 ──
    phase_data = {}
    phase_ct = {}
    for p in phase_keys:
        seg_img = nib.load(phases[p]["seg"])
        phase_data[p] = np.round(np.asanyarray(seg_img.dataobj)).astype(np.uint16)
        img_path = phases[p]["img"]
        if img_path and os.path.exists(img_path):
            phase_ct[p] = np.asanyarray(nib.load(img_path).dataobj).astype(np.float32)
        else:
            phase_ct[p] = None

    label_names = {0: "배경", 1: "신장", 2: "종양", 3: "물혹"}

    # Shape 정보 출력
    print(f"\n  Phase별 Shape:")
    for p in phase_keys:
        print(f"    {p}: {phase_data[p].shape}")

    # ── 2. Phase별 라벨 크기/비율 ──
    print(f"\n{'='*60}")
    print(f"  Phase별 라벨 크기 비교")
    print(f"{'='*60}")

    # 헤더
    header = f"  {'라벨':<10}"
    for p in phase_keys:
        header += f"  {p:>12}"
    header += f"  {'최대차이':>10}"
    print(header)
    print(f"  {'─'*10}" + f"  {'─'*12}" * len(phase_keys) + f"  {'─'*10}")

    for label in [1, 2, 3]:
        name = label_names[label]
        counts = []
        for p in phase_keys:
            counts.append(int(np.sum(phase_data[p] == label)))

        if all(c == 0 for c in counts):
            continue

        row = f"  {name:<10}"
        for c in counts:
            row += f"  {c:>12,}"

        max_c = max(counts) if max(counts) > 0 else 1
        min_c = min(counts)
        diff_pct = (max_c - min_c) / max_c * 100
        row += f"  {diff_pct:>9.1f}%"
        print(row)

    # 전체 장기
    row = f"  {'장기전체':<10}"
    organ_counts = []
    for p in phase_keys:
        c = int(np.sum(phase_data[p] > 0))
        organ_counts.append(c)
        row += f"  {c:>12,}"
    max_c = max(organ_counts) if max(organ_counts) > 0 else 1
    min_c = min(organ_counts)
    diff_pct = (max_c - min_c) / max_c * 100
    row += f"  {diff_pct:>9.1f}%"
    print(row)

    # ── 3. Phase 쌍별 Dice 계수 ──
    print(f"\n{'='*60}")
    print(f"  Phase 쌍별 Dice 계수")
    print(f"{'='*60}")

    from itertools import combinations
    for p1, p2 in combinations(phase_keys, 2):
        print(f"\n  [{p1} vs {p2}]")

        if phase_data[p1].shape != phase_data[p2].shape:
            print(f"    ⚠ Shape 불일치: {p1}={phase_data[p1].shape}, {p2}={phase_data[p2].shape} — 건너뜀")
            continue

        for label in [1, 2, 3]:
            name = label_names[label]
            mask1 = (phase_data[p1] == label)
            mask2 = (phase_data[p2] == label)
            sum1 = int(np.sum(mask1))
            sum2 = int(np.sum(mask2))

            if sum1 == 0 and sum2 == 0:
                continue

            intersection = int(np.sum(mask1 & mask2))
            dice = 2 * intersection / (sum1 + sum2) if (sum1 + sum2) > 0 else 0.0
            print(f"    {name}: Dice={dice:.4f}  (겹침={intersection:,}, {p1}={sum1:,}, {p2}={sum2:,})")

        # 장기 전체
        mask1 = (phase_data[p1] > 0)
        mask2 = (phase_data[p2] > 0)
        intersection = int(np.sum(mask1 & mask2))
        sum1 = int(np.sum(mask1))
        sum2 = int(np.sum(mask2))
        dice = 2 * intersection / (sum1 + sum2) if (sum1 + sum2) > 0 else 0.0
        print(f"    장기전체: Dice={dice:.4f}  (겹침={intersection:,}, {p1}={sum1:,}, {p2}={sum2:,})")

    # ── 4. 불일치 영역 슬라이스 분석 ──
    print(f"\n{'='*60}")
    print(f"  불일치 영역 슬라이스 분석 (axial 기준)")
    print(f"{'='*60}")

    for p1, p2 in combinations(phase_keys, 2):
        print(f"\n  [{p1} vs {p2}]")

        if phase_data[p1].shape != phase_data[p2].shape:
            print(f"    ⚠ Shape 불일치 — 건너뜀")
            continue

        for label in [1, 2, 3]:
            name = label_names[label]
            mask1 = (phase_data[p1] == label)
            mask2 = (phase_data[p2] == label)

            if int(np.sum(mask1)) == 0 and int(np.sum(mask2)) == 0:
                continue

            # 슬라이스별 차이 계산 (axial = 마지막 축 기준이 아닌 첫 축 기준)
            n_slices = mask1.shape[2]
            diff_per_slice = np.zeros(n_slices)
            for s in range(n_slices):
                slice1 = mask1[:, :, s]
                slice2 = mask2[:, :, s]
                diff_per_slice[s] = int(np.sum(slice1 != slice2))

            total_diff = int(np.sum(diff_per_slice > 0))
            if total_diff == 0:
                print(f"    {name}: 완전 일치")
                continue

            # 차이가 있는 슬라이스 범위
            diff_slices = np.where(diff_per_slice > 0)[0]
            # 상위 차이 슬라이스
            top_slices = np.argsort(diff_per_slice)[::-1][:5]
            top_slices = top_slices[diff_per_slice[top_slices] > 0]

            print(f"    {name}: 불일치 슬라이스 {total_diff}개 (범위: {diff_slices[0]}~{diff_slices[-1]})")
            print(f"      가장 큰 차이:")
            for s in top_slices:
                only1 = int(np.sum(mask1[:, :, s] & ~mask2[:, :, s]))
                only2 = int(np.sum(mask2[:, :, s] & ~mask1[:, :, s]))
                print(f"        slice {s}: {p1}에만 {only1:,}, {p2}에만 {only2:,} voxels")

    print()


# ──────────────────────────────────────────────
# 영역 지정 Region Restriction
# ──────────────────────────────────────────────

def build_region_mask(shape, affine):
    """슬라이스 범위, 바운딩 박스, 방향을 조합하여 3D 마스크 생성."""
    mask = np.ones(shape, dtype=bool)

    print("\n  영역 지정 모드 Region Restriction Mode")
    print("  현재 shape Current shape:", shape)

    # ── 1. 슬라이스 범위 Slice range ──
    use_slice = input_choice("  슬라이스 범위 제한? Restrict slice range?", [
        "1: 예 Yes",
        "2: 아니오 No (전체 슬라이스 All slices)",
    ])
    if use_slice == "1":
        axis = input_choice("  축 선택 Select axis", [
            "0: axis 0 (shape={})".format(shape[0]),
            "1: axis 1 (shape={})".format(shape[1]),
            "2: axis 2 (shape={})".format(shape[2]),
        ])
        axis = int(axis)
        start = input_int(f"  시작 슬라이스 Start slice (기본 default 0)", default=0)
        end = input_int(f"  끝 슬라이스 End slice (기본 default {shape[axis]-1})", default=shape[axis]-1)
        start = max(0, min(start, shape[axis]-1))
        end = max(start, min(end, shape[axis]-1))
        slicing = [slice(None)] * 3
        slicing[axis] = slice(0, start)
        mask[tuple(slicing)] = False
        slicing[axis] = slice(end + 1, shape[axis])
        mask[tuple(slicing)] = False
        print(f"  → axis {axis}, 슬라이스 slice {start}~{end} ({end-start+1} slices)")

    # ── 2. 3D 바운딩 박스 Bounding box ──
    use_bbox = input_choice("  바운딩 박스 제한? Restrict bounding box?", [
        "1: 예 Yes",
        "2: 아니오 No",
    ])
    if use_bbox == "1":
        for ax, name in enumerate(["X (axis 0)", "Y (axis 1)", "Z (axis 2)"]):
            lo = input_int(f"  {name} 시작 start (기본 default 0)", default=0)
            hi = input_int(f"  {name} 끝 end (기본 default {shape[ax]-1})", default=shape[ax]-1)
            lo = max(0, min(lo, shape[ax]-1))
            hi = max(lo, min(hi, shape[ax]-1))
            slicing = [slice(None)] * 3
            slicing[ax] = slice(0, lo)
            mask[tuple(slicing)] = False
            slicing[ax] = slice(hi + 1, shape[ax])
            mask[tuple(slicing)] = False
            print(f"  → {name}: {lo}~{hi}")

    # ── 3. 방향 제한 Direction restriction ──
    use_dir = input_choice("  방향 제한? Restrict direction?", [
        "1: 예 Yes",
        "2: 아니오 No",
    ])
    if use_dir == "1":
        ax_codes = nib.aff2axcodes(affine)
        print(f"  축 방향 Axis orientation: axis0={ax_codes[0]}, axis1={ax_codes[1]}, axis2={ax_codes[2]}")

        # 각 축의 양쪽 방향 매핑
        opposites = {'R': 'L', 'L': 'R', 'A': 'P', 'P': 'A', 'S': 'I', 'I': 'S'}
        dir_names = {
            'R': 'Right 우', 'L': 'Left 좌',
            'A': 'Anterior 전', 'P': 'Posterior 후',
            'S': 'Superior 상', 'I': 'Inferior 하',
        }

        # 6방향 중 선택
        all_dirs = []
        for ax_idx, code in enumerate(ax_codes):
            opp = opposites[code]
            all_dirs.append((code, ax_idx, "low"))
            all_dirs.append((opp, ax_idx, "high"))

        options = [f"{i+1}: {dir_names[d[0]]} ({d[0]}, axis {d[1]})" for i, d in enumerate(all_dirs)]
        dir_choice = input_choice("  제한할 방향 Direction to restrict", options)
        chosen = all_dirs[int(dir_choice) - 1]
        dir_code, dir_axis, dir_side = chosen

        mid = shape[dir_axis] // 2
        cut = input_int(
            f"  기준 슬라이스 Cut at slice (기본 default {mid}, 범위 range 0~{shape[dir_axis]-1})",
            default=mid)
        cut = max(0, min(cut, shape[dir_axis]-1))

        slicing = [slice(None)] * 3
        if dir_side == "low":
            # 해당 방향 = 낮은 인덱스 쪽 → 낮은 쪽만 유지
            slicing[dir_axis] = slice(cut, shape[dir_axis])
            mask[tuple(slicing)] = False
            print(f"  → {dir_names[dir_code]} 방향: axis {dir_axis}, slice 0~{cut-1} 유지")
        else:
            # 해당 방향 = 높은 인덱스 쪽 → 높은 쪽만 유지
            slicing[dir_axis] = slice(0, cut + 1)
            mask[tuple(slicing)] = False
            print(f"  → {dir_names[dir_code]} 방향: axis {dir_axis}, slice {cut+1}~{shape[dir_axis]-1} 유지")

    restricted = int(np.sum(mask))
    total = int(np.prod(shape))
    print(f"\n  영역 지정 완료 Region set: {restricted:,} / {total:,} voxels ({restricted/total*100:.1f}%)")
    return mask


def apply_with_region(func, data, region_mask, **kwargs):
    """영역 지정 마스크 적용: 전체 데이터로 실행 후 영역 안 변경만 반영."""
    result = func(data, **kwargs)
    # 변경된 복셀 중 영역 안 + 영역 경계 1복셀 확장 범위만 반영
    struct = ndimage.generate_binary_structure(3, 1)
    expanded_region = ndimage.binary_dilation(region_mask, structure=struct, iterations=1)
    final = data.copy()
    final[expanded_region] = result[expanded_region]
    return final


# ──────────────────────────────────────────────
# 메인 루프
# ──────────────────────────────────────────────

FUNCTIONS = {
    "1": ("라벨 상태 분석 Label Analysis", func_analyze),
    "2": ("고립 복셀 제거 Remove Isolated", func_remove_isolated),
    "3": ("저강도 제거 Remove Low Intensity (≤ 0)", func_remove_low_intensity),
    "4": ("고강도 제거 Remove High Intensity (≥ threshold)", func_remove_high_intensity),
    "5": ("Smoothing (종양 Tumor/물혹 Cyst/장기 전체 Organ)", func_smooth),
    "6": ("경계 확장 Boundary Expansion (신장/종양/물혹)", func_expand),
    "7": ("경계 축소 Boundary Trimming (장기 외곽 Organ/종양/물혹)", func_trim_boundary),
    "8": ("경계 계단 메꿈 Staircase Fill", func_fill_staircase),
    "9": ("돌출부 제거 Protrusion Removal", func_remove_protrusion),
    "10": ("내부 구멍 채우기 Fill Internal Holes", func_fill_holes),
    "11": ("고립 신장→종양 재라벨링 Relabel Isolated Kidney→Tumor", func_relabel_isolated_kidney),
    "12": ("볼록 껍질 라벨링 Convex Hull Labeling (종양/물혹)", func_label_convex),
    "13": ("세그멘테이션 합치기 Merge Segmentations", func_merge_segmentations),
    "14": ("Phase 비교 분석 Phase Comparison", None),
    "m": ("영역 지정 실행 Region Restricted Run", None),
    "r": ("롤백 Rollback", None),
}


def main():
    if len(sys.argv) < 2:
        print("사용법: python segtools.py <케이스 폴더>")
        print("예시:   python segtools.py S004")
        sys.exit(1)

    case_dir = sys.argv[1]
    if not os.path.isdir(case_dir):
        print(f"폴더를 찾을 수 없음: {case_dir}")
        sys.exit(1)

    case_name = os.path.basename(os.path.abspath(case_dir))
    phases = load_case(case_dir)

    if not phases:
        print(f"세그멘테이션 파일을 찾을 수 없음: {case_dir}")
        sys.exit(1)

    print(f"\n케이스 Case: {case_name}")
    print(f"사용 가능한 phase Available phases: {', '.join(sorted(phases.keys()))}")

    while True:
        # Phase 선택
        print(f"\n{'─'*50}")
        print(f"Phase 선택 Select phase (q: 종료 quit)")
        for p in sorted(phases.keys()):
            ct_status = "CT있음" if phases[p]["img"] else "CT없음"
            print(f"  {p}: {os.path.basename(phases[p]['seg'])} ({ct_status})")
        print(f"  all: 모든 phase에 동일 작업 실행 Apply to all phases")

        phase_input = input("\n  Phase: ").strip()
        if phase_input.upper() in ("Q", "QUIT", "EXIT", "ㅂ"):
            print("종료.")
            break

        # 한글 입력 매핑
        hangul_phase = {"ㅁ": "A", "ㅇ": "D", "ㅔ": "P"}
        if phase_input in hangul_phase:
            phase_input = hangul_phase[phase_input]

        if phase_input.lower() == "all":
            selected_phases = sorted(phases.keys())
        elif phase_input.upper() in phases:
            selected_phases = [phase_input.upper()]
        else:
            print("  → 유효한 phase를 입력하세요")
            continue

        # 기능 선택
        print(f"\n  기능 선택 Select function:")
        for key, (name, _) in FUNCTIONS.items():
            # 롤백은 히스토리 있을 때만 표시
            if key == "r":
                has_history = any(p in rollback_history and len(rollback_history[p]) > 0
                                 for p in selected_phases)
                if has_history:
                    print(f"    {key}: {name}")
            else:
                print(f"    {key}: {name}")
        print(f"    b: 돌아가기 Back")

        func_input = input("\n  기능 Function: ").strip().lower()
        if func_input in ("b", "ㅠ"):
            continue
        if func_input not in FUNCTIONS:
            print("  → 유효한 기능 번호를 입력하세요")
            continue

        func_name, func = FUNCTIONS[func_input]

        # ── Phase 비교 분석 처리 ──
        if func_input == "14":
            func_compare_phases(phases)
            continue

        # ── 영역 지정 실행 처리 ──
        if func_input in ("m", "ㅡ"):
            # 영역 지정 가능한 기능 목록 (분석/비교/롤백/영역지정 자체 제외)
            region_funcs = {k: v for k, v in FUNCTIONS.items()
                           if k not in ("1", "14", "m", "r") and v[1] is not None}

            print(f"\n  영역 지정 후 실행할 기능 선택 Select function for region run:")
            for key, (name, _) in region_funcs.items():
                print(f"    {key}: {name}")
            print(f"    b: 돌아가기 Back")

            rf_input = input("\n  기능 Function: ").strip().lower()
            if rf_input in ("b", "ㅠ") or rf_input not in region_funcs:
                if rf_input != "b":
                    print("  → 유효한 기능 번호를 입력하세요")
                continue

            rf_name, rf_func = region_funcs[rf_input]

            for phase in selected_phases:
                seg_path = phases[phase]["seg"]
                img_path = phases[phase]["img"]

                print(f"\n{'='*50}")
                print(f"  [{phase} phase] 영역 지정 Region Restriction + {rf_name}")
                print(f"{'='*50}")

                seg_img = nib.load(seg_path)
                data = np.round(np.asanyarray(seg_img.dataobj)).astype(np.uint16)
                zooms = seg_img.header.get_zooms()

                ct_data = None
                if img_path and os.path.exists(img_path):
                    ct_img = nib.load(img_path)
                    ct_data = np.asanyarray(ct_img.dataobj).astype(np.float32)

                # 영역 지정
                try:
                    region_mask = build_region_mask(data.shape, seg_img.affine)
                except CancelOperation:
                    print("\n  취소됨 Cancelled.")
                    continue

                # 기능 실행 (영역 제한)
                try:
                    result = apply_with_region(rf_func, data, region_mask,
                                               ct_data=ct_data, zooms=zooms)
                except CancelOperation:
                    print("\n  취소됨 Cancelled.")
                    continue

                if np.array_equal(data, result):
                    print("\n  변경 없음 No changes.")
                    continue

                changed = int(np.sum(data != result))
                print(f"\n  변경된 복셀 Changed voxels: {changed:,}")

                if phase not in rollback_history:
                    rollback_history[phase] = []
                rollback_history[phase].append((data.copy(), f"[영역지정] {rf_name}", seg_img))

                backup_file(seg_path)
                save_result(seg_path, result, seg_img)

                history_depth = len(rollback_history[phase])
                print(f"  (롤백 가능 Rollback available: {history_depth}단계 steps)")

                while True:
                    prompt = "  다음 Next? (enter:계속 continue / r:롤백 rollback"
                    if len(selected_phases) > 1:
                        prompt += " / q:나머지 건너뛰기 skip rest"
                    prompt += "): "
                    next_input = input(prompt).strip().lower()
                    if next_input == "":
                        break
                    elif next_input in ("q", "ㅂ") and len(selected_phases) > 1:
                        break
                    elif next_input in ("r", "ㄱ"):
                        if phase in rollback_history and len(rollback_history[phase]) > 0:
                            prev_data, prev_desc, prev_img = rollback_history[phase][-1]
                            print(f"  [{phase} phase] 롤백: '{prev_desc}' 이전 상태로 되돌리기")
                            save_result(seg_path, prev_data, prev_img)
                            rollback_history[phase].pop()
                            print(f"  남은 히스토리: {len(rollback_history[phase])}단계")
                        else:
                            print(f"  [{phase} phase] 롤백 히스토리 없음")
                        break
                    else:
                        valid = "enter / r / q" if len(selected_phases) > 1 else "enter / r"
                        print(f"    → {valid} 중 하나를 입력하세요")
                if next_input in ("q", "ㅂ"):
                    break
            continue

        # ── 롤백 처리 ──
        if func_input in ("r", "ㄱ"):
            for phase in selected_phases:
                seg_path = phases[phase]["seg"]
                if phase not in rollback_history or len(rollback_history[phase]) == 0:
                    print(f"\n  [{phase} phase] 롤백 히스토리 없음")
                    continue

                prev_data, prev_desc, seg_img = rollback_history[phase][-1]
                print(f"\n  [{phase} phase] 롤백: '{prev_desc}' 이전 상태로 되돌리기")

                save_result(seg_path, prev_data, seg_img)
                rollback_history[phase].pop()

                remain = len(rollback_history[phase])
                print(f"  남은 히스토리: {remain}단계")
            continue

        # ── 일반 기능 실행 ──
        for phase in selected_phases:
            seg_path = phases[phase]["seg"]
            img_path = phases[phase]["img"]

            print(f"\n{'='*50}")
            print(f"  [{phase} phase] {func_name}")
            print(f"{'='*50}")

            # 로드
            seg_img = nib.load(seg_path)
            data = np.round(np.asanyarray(seg_img.dataobj)).astype(np.uint16)
            zooms = seg_img.header.get_zooms()

            ct_data = None
            if img_path and os.path.exists(img_path):
                ct_img = nib.load(img_path)
                ct_data = np.asanyarray(ct_img.dataobj).astype(np.float32)

            # 실행
            try:
                result = func(data, ct_data=ct_data, zooms=zooms)
            except CancelOperation:
                print("\n  취소됨 Cancelled.")
                continue

            # 변경 확인
            if np.array_equal(data, result):
                print("\n  변경 없음 No changes.")
                continue

            changed = int(np.sum(data != result))
            print(f"\n  변경된 복셀 Changed voxels: {changed:,}")

            # 저장 (히스토리에 이전 상태 저장 후 결과 저장)
            if phase not in rollback_history:
                rollback_history[phase] = []
            rollback_history[phase].append((data.copy(), func_name, seg_img))

            backup_file(seg_path)
            save_result(seg_path, result, seg_img)

            history_depth = len(rollback_history[phase])
            print(f"  (롤백 가능 Rollback available: {history_depth}단계 steps)")

            # 저장 후 확인
            while True:
                prompt = "  다음 Next? (enter:계속 continue / r:롤백 rollback"
                if len(selected_phases) > 1:
                    prompt += " / q:나머지 건너뛰기 skip rest"
                prompt += "): "
                next_input = input(prompt).strip().lower()
                if next_input == "":
                    break
                elif next_input in ("q", "ㅂ") and len(selected_phases) > 1:
                    break
                elif next_input in ("r", "ㄱ"):
                    if phase in rollback_history and len(rollback_history[phase]) > 0:
                        prev_data, prev_desc, prev_img = rollback_history[phase][-1]
                        print(f"  [{phase} phase] 롤백: '{prev_desc}' 이전 상태로 되돌리기")
                        save_result(seg_path, prev_data, prev_img)
                        rollback_history[phase].pop()
                        print(f"  남은 히스토리: {len(rollback_history[phase])}단계")
                    else:
                        print(f"  [{phase} phase] 롤백 히스토리 없음")
                    break
                else:
                    valid = "enter / r / q" if len(selected_phases) > 1 else "enter / r"
                    print(f"    → {valid} 중 하나를 입력하세요")
            if next_input in ("q", "ㅂ"):
                break


if __name__ == "__main__":
    main()
