# Performance Benchmarks

Cross-version performance measurements for Ptera Software's solvers. Two probes are charted: solver wall time and solver peak memory delta. The data is produced by a separate benchmark archive that runs each released version against the same configuration suite on a fixed host. That archive is kept private because its self-hosted runner introduces potential security concerns. If you are interested in working on Ptera Software's performance and would like access, email <camerongurban@gmail.com>.

## Methodology

Each version's bar is measured against a frozen software stack reconstructed at that version's release day. Concretely, the Python interpreter (exact patch) and the dependency set (pinned with content hashes) are the ones that a fresh install on the release day would have resolved. The bars therefore reflect release-day resolution held constant, not present-day resolution. Installs use `uv` with `--require-hashes` so the install path is byte-identical across re-runs.

Every bar in the archive shares the host hardware described under [Test environment](#test-environment) below. The publishing step refuses to emit a chart that mixes bars from different hosts, so a hardware change in the archive forces a full re-run before publication can resume.

The baseline against which speedups and ratios are computed is the immediately-preceding pypi version of the latest pypi version. Open the interactive chart's toggle ("Speedup vs. baseline" for time, "Ratio vs. baseline" for memory) to switch from absolute units on a log y-axis to the normalized view on a linear y-axis.

## Test environment

```{list-table}
:header-rows: 0
:widths: 20 80

* - OS
  - {{ host_os }}
* - CPU
  - {{ host_cpu }} ({{ host_cores }} logical cores, governor: {{ host_governor }})
* - Memory
  - {{ host_memory_mb }} MB (swappiness: {{ host_swappiness }}, transparent hugepages: {{ host_thp }})
* - GPU
  - {{ host_gpu }} (driver {{ host_gpu_driver }}, CUDA {{ host_cuda }})
```

The full redacted host snapshot is also published as <a href="_static/benchmarks/host.json">host.json</a>.

## Current results

```{raw} html
<iframe title="Solver run time across versions (light theme)"
        src="_static/benchmarks/solver_run_time_light.html"
        class="only-light"
        width="100%"
        height="700"
        style="border: 0; display: block;"
        loading="lazy"></iframe>
<iframe title="Solver run time across versions (dark theme)"
        src="_static/benchmarks/solver_run_time_dark.html"
        class="only-dark"
        width="100%"
        height="700"
        style="border: 0; display: block;"
        loading="lazy"></iframe>
```

Download the raw data: <a href="_static/benchmarks/solver_run_time.csv">solver_run_time.csv</a>

```{raw} html
<iframe title="Solver run peak memory delta across versions (light theme)"
        src="_static/benchmarks/solver_run_memory_light.html"
        class="only-light"
        width="100%"
        height="700"
        style="border: 0; display: block;"
        loading="lazy"></iframe>
<iframe title="Solver run memory delta across versions (dark theme)"
        src="_static/benchmarks/solver_run_memory_dark.html"
        class="only-dark"
        width="100%"
        height="700"
        style="border: 0; display: block;"
        loading="lazy"></iframe>
```

Download the raw data: <a href="_static/benchmarks/solver_run_memory.csv">solver_run_memory.csv</a>

## Notes on what's shown

- Configurations are ordered along the x-axis by cost. The slug prefix (for example, `0007_8E07evals_urvlm`) records the approximate Biot-Savart kernel-evaluation count and the solver family (`srvlm` for the steady ring vortex lattice method, `shvlm` for the steady horseshoe vortex lattice method, or `urvlm` for the unsteady ring vortex lattice method).
- Cells deliberately gated out of a version (for example, configurations too large to fit in memory on this host) render as missing bars rather than as zero-height bars. The CSVs preserve the distinction between `ok`, `skipped`, and `missing` states in the `status` column, with the reason (when known) in the `reason` column.
- No aggregate or "headline" speedup is reported. The benchmark suite is a trimmed, down-sampled selection from the Monte Carlo distribution of expected workloads, so an arithmetic or geometric mean over the selection would not faithfully represent the underlying distribution. The chart is the report; read individual configurations.
