# Validation scripts

The Python scripts in this directory back the CTest release suite.

Recommended interface:

```bash
ctest --test-dir build --output-on-failure
```

Individual checks can also be run directly, for example:

```bash
python validation/run_smoke.py --exe build/nanokmc
python validation/run_output_smoke.py --exe build/nanokmc
python validation/run_profile_axis_smoke.py --exe build/nanokmc
python validation/run_rate_category_limits.py --exe build/nanokmc
python validation/run_input_validation.py --exe build/nanokmc
```

See [`../docs/validation.md`](../docs/validation.md) for the scientific invariants and scope of each check.
