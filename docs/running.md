# Running the pipeline

Run from the repository root with the controller environment active. Commands
below use a complete project saved as `config/data/local.yaml`.

## Inspect, then execute locally

```bash
snakemake --configfile config/data/local.yaml \
  --cores 8 --use-conda --conda-frontend conda --dry-run --printshellcmds
snakemake --configfile config/data/local.yaml \
  --cores 8 --use-conda --conda-frontend conda --printshellcmds --rerun-incomplete
```

Initialization creates `out_dir` and writes logs, processed configurations, and
`requested_files.tab`, including during dry runs. Set `out_dir` and `tmp_dir`
together; only one workflow should write to a given output tree at a time.

## Target a file or extend a stage chain

Stage names such as `contig` are configuration entries. For an individual result,
use its path from `requested_files.tab`, for example:

```bash
snakemake --configfile config/data/local.yaml \
  --cores 8 --use-conda --dry-run \
  results/sample/contig/hifiasm_l1_no_ul@p2/sample.contig.hap1.len
```

This path assumes `out_dir: results/sample/`, `genome_prefix: sample`, and the
branch `hifiasm_l1_no_ul@p2`. To extend an assembly, retain upstream stages and
presets in `stage_list`; see [stage ordering](stages.md#build-an-assembly-chain-incrementally).

## Resources

| Pipeline setting | Effect |
| --- | --- |
| `threads.tab` / rule `threads` | Threads supplied to tool commands; some rules also use separate thread parameters. |
| `memory_mb.tab` / `resources.mem` | Memory estimate in MB, forwarded as `sbatch --mem` by the legacy Slurm profile. |
| `time.tab` / `resources.time` | Scheduler time request. |
| `resources.cpus` | CPU allocation forwarded by the Slurm submission string; separate from rule `threads`. |

The `get_threads` helper caps `resources.cpus` only in rules that call it; it does
not replace rule `threads`. Check both tool threads and scheduler CPU requests.
Neither `--cores` nor the configured queue `max_mem` imposes a local RAM limit.
Some rules increase memory requests on retries; the legacy profile sets
`retries: 0`.

For persistent changes, use private copies of all three resource tables; see
[resource configuration](configuration.md#resource-profiles-and-parameter-sets).

## Slurm: adapt the legacy profile

For Snakemake 7, use `profile/slurm/`. Cluster settings are loaded from
`config/cluster.yaml` through `cluster_config_file` in the core configuration.
Adapt these settings to your cluster:

- `queue.cpu.name`, queue limits, `nodes`, and per-rule `rule_nodes`.
- The profile's `cluster` submission string, including `--qos=normal` and site
  account options. It maps `queue`, `cpus`, `mem`, and `time` to Slurm requests.
- Shared paths and tool environments accessible to both compute nodes and the
  controller; rules marked `localrules` execute on the controller.

```bash
snakemake --profile profile/slurm \
  --configfile config/data/local.yaml --jobs 20 --local-cores 2 \
  --conda-frontend conda --dry-run --printshellcmds
```

The profile enables Conda and incomplete-job reruns. The alternative
`profile/slurm8.30/` uses an executor-style configuration and is not a verified
replacement for this Snakemake 7 profile. See the
[Snakemake CLI reference](https://snakemake.readthedocs.io/en/v7.32.0/executing/cli.html)
for generic execution and recovery options.
