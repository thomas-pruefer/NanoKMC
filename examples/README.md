# Examples

The repository ships one small run directory for each public solver:

- `classical` — `KMCClassical`
- `active_filtered_generic` — `KMCActiveFilteredGeneric`
- `active_filtered_binary_nn` — `KMCActiveFilteredBinaryNN`
- `partial_filter_optimized` — `KMCPartialFilterOptimized`
- `rate_category_optimized` — `KMCRateCategoryOptimized`
- `exact_class_optimized` — `KMCExactClassOptimized`

An additional `output_features` example demonstrates the documented output surface using the Active-Filtered BinaryNN solver.

Each run directory contains `nanokmc.in` and the required `output/CalcData.csv`. NanoKMC creates evaluation and checkpoint subdirectories as needed. `CalcData.template.csv` provides a clean checkpoint file for resetting an example after a run.
