
# Adamselv Hydrodynamic Model Setup – Documentation & User Guide

## Project Structure Overview

```
Adamselv/
├── 01_build_grid/
├── 02_build_rivers/
├── 03_atm_forcing/
├── 04_add_fabm_tracer/
├── 05_qc/
```

---

## 1. Grid Generation (`01_build_grid/`)

### Purpose:
This step builds the model grid from 2DM and bathymetric data, preparing the simulation domain.

### Key Scripts:
- `2dm_to_grid_fvtools.py`: Converts SMS 2DM files into structured grid data.
- `boundary_prepro.py`: Processes boundary conditions, like open ocean boundaries.
- `make_grd_adamselv.py`: Main driver script for assembling the grid and domain configuration.

### Inputs:
- `.2dm` files (e.g., `adamselv_v01.2dm`)
- Bathymetry data (`bathy_hybrid_adamselv_v01.npy`)
- Depth file (`adamselv_dep_NEW.dat`)

### Outputs:
- Grid files and initial configuration `.nc` and `.npy` files.

> ⚠️ **Important**: After the grid is built, a restart file must be generated on the server. This restart is crucial to initialize the model with the built grid and prepare it for simulation.

---

## 2. River Forcing (`02_build_rivers/`)

### Purpose:
Add river inflows and process their grid representation.

### Key Scripts:
- `build_rivers.py`: Builds river inflow boundary files.
- `convert_to_node.py`: Converts river input locations to grid nodes.

### Inputs:
- River discharge/temperature data (e.g., `Finnmark_temperatures.npy`)

### Outputs:
- NetCDF or text-based boundary condition inputs for rivers.

---

## 3. Atmospheric Forcing (`03_atm_forcing/`)

### Purpose:
Processes and fixes atmospheric data (e.g., evaporation and precipitation) used to force the model.

### Key Scripts:
- `fix_evaporation.py`: Cleans and adjusts evaporation data.
- `check_MEPS_evap.py`: Quality checks MEPS model evaporation outputs.
- `run_effie.py`: Wrapper or main driver to process meteorological forcing.

### Outputs:
- Processed atmospheric data in `output/`
- Logs such as `log_effie.txt`

---

## 4. FABM Tracer & Restart Files (`04_add_fabm_tracer/`)

### Purpose:
This step integrates the FABM biogeochemical model and tracer data into the hydrodynamic model.

### Key Components:
- `add_fabm_river.ipynb`: Jupyter notebook to add FABM-related fields and inputs.
- `make_input.py`, `make_all_water_col.py`: Python scripts to build water column inputs.

### Outputs:
- Restart and tracer-ready NetCDF files
- FABM-compatible river files

---

## 5. Quality Control (`05_qc/`)

### Purpose:
Run quality checks on grid and atmospheric data to ensure consistency before simulation.

### Scripts:
- `qc_output.py`: Validates output files.
- `qc_atm_forcing.py`: Checks integrity of atmospheric data.

---

## Running the Model (Simplified Workflow)

1. **Build the Grid**
   ```bash
   python make_grd_adamselv.py
   ```

2. **Build Rivers**
   ```bash
   python build_rivers.py
   python convert_to_node.py
   ```

3. **Process Atmospheric Forcing**
   ```bash
   python fix_evaporation.py
   python check_MEPS_evap.py
   python run_effie.py
   ```

4. **⚠️ Generate Restart on Server**
   - This step is not fully scripted here and must be done after grid generation.
   - Upload or move generated grid to the server.
   - Use model engine (e.g., FVCOM or similar) to generate a restart file based on the initial grid.

5. **Add FABM Tracers**
   - Use Jupyter notebook `add_fabm_river.ipynb` to configure and add tracers.
   - Run supporting scripts to finalize input.

6. **Run Quality Checks**
   ```bash
   python qc_output.py
   python qc_atm_forcing.py
   ```

---

## Notes

- Make sure dependencies such as NumPy, netCDF4, and any custom grid libraries (like FVTools) are installed.
- Always verify file paths and environment variables, especially when running on a cluster or server.
