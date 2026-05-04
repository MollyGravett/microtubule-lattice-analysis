# microtubule-lattice-analysis
## Requirements

- [IMOD](https://bio3d.colorado.edu/imod/)
- [PEET](https://bio3d.colorado.edu/PEET/)
- [EMAN2](https://cryoem.bcm.edu/cryoem/downloads/view_eman2_versions)
- [RELION](https://relion.readthedocs.io/)

## Microtubule Protofilament References

References available from:
[https://github.com/moores-lab/MiRPv2/tree/main/data/protofilament_sorting_references](https://github.com/moores-lab/MiRPv2/tree/main/data/protofilament_sorting_references)

---

## Pipeline Steps


### 1. Reconstruct Tomograms

Reconstruct tomograms in AreTomo or IMOD.

---

### 2. Pick Coordinates in IMOD

Pick coordinates in IMOD and save as a `.mod` file (each microtubule as a new contour).

---

### 3. Unbinned Aligned Stack

```bash
subm newst.com
```

---

### 4. Unbinned Tomogram Reconstruction

```bash
mousepad tilt.com
```

Edit `tilt.com`:
- Comment out `LOCALFILE` and `XTILTFILE`
- Replace `PERPENDICULAR` with `RotateBy90`

```bash
subm tilt.com
```

---

### 5. Create Unbinned `.mod`

```bash
imodtrans -2 AreTomo3Output/WT_tomo01.mrc_Imod/WT_tomo01.mrc_st.xf \
  -I AreTomo3Output/WT_tomo01.mrc_Imod/WT_tomo01.mrc_st.ali \
  WT_tomo01.mod disco_tomo01_unbinned.mod
```

---

### 6. PEET Add Mod Points

```bash
addModPts WT_tomo01_unbinned.mod 9.5
# Dimer distance in pixels, adjust depending on pixel size
```

---

### 7. Open with Unbinned Tomo and Resave

```bash
3dmod -m WT_tomo01.mrc_st.rec WT_tomo01_unbinned_PtsAdded.mod
# Save as: WT_tomo01_unbinned_PtsAdded_resaved.mod
```

---

### 8. Extract Subtomograms

```bash
mkdir subtomos_unbinned

boxstartend \
  -image disco_tomo01.mrc_st.rec \
  -model WT_tomo01_unbinned_resaved_PtsAdded.mod \
  -series subtomos_unbinned/WT_tomo01_subtomo \
  -box 300
```
---

### 9. Rescale MiRP References to Match Subtomograms

```bash
relion_helix_toolbox \
  --i 11pf_syn_ref_tubulin_only_6A_5-56Apix.mrc \
  --o 11pf_syn_ref_tubulin_only_lpf20_2p21Apix_box_300.mrc \
  --rescale_angpix 2.21 \
  --new_box 300 \
  --lowpass 20
```

---

### 10. Make Tube with RELION Helical Toolbox

```bash
relion_helix_toolbox --cylinder \
  --o ptcl_sub_tube_ref_box300_2p21_out300.mrc \
  --boxdim 300 \
  --cyl_outer_diameter 300 \
  --angpix 2.21
```

---

### 11. Make Soft Mask of Tube in RELION

Run a **Mask Create** job in RELION.

---

### 12. Invert References

Edit `invert.py`, then run:

```bash
python invert.py
```

---

### 13. Make Reference Stack with Rotational Averages

> Requires EMAN2

Edit `reference_stack_rot_av.py`, then run:

```bash
python reference_stack_rot_av.py
```

---

### 14. Make RELION 3 Particle `.star` from `.mod`

```bash
python imod_model_to_rln3_star_file_per_MT.py
```

---

### 15. Measure Lattice Spacing

```bash
mkdir lattice_spacing
cd lattice_spacing
```

> Requires EMAN2

```bash
python padded_layer_line_big_crop_pf_polarity.py
```

Measures layer lines, orientation, and protofilament number.

