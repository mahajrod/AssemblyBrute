# Configuration

Start with a complete project configuration from the
[example catalog](examples/README.md) and adapt it to your sample.
The fragments below illustrate individual settings; they are not standalone
project files. Run commands from the repository root.

## Which file to edit

| File or setting | Purpose |
| --- | --- |
| Project YAML passed to `--configfile` | Sample identity, datatypes, stages, tool/preset selections, QC levels, and explicit skip flags. |
| Core YAML selected by `main_config_file` | Shared configuration, path settings, environment mappings, configuration-file locations, and parameter defaults/overrides. |
| `config/coretools/<parameter_set>/<tool>.tab` | Named Hifiasm, Flye, and Verkko presets. Shipped tables are under `default/`. |
| `config/core.secondary_tools.yaml` | Default selections for secondary-tool option sets. |
| `config/resources/<resource_profile>/` | `threads.tab`, `memory_mb.tab`, and `time.tab`, indexed by rule and parameter-set column. |
| `config/cluster.yaml` | Queue names, queue limits, and node selections. This is pipeline configuration, not a Snakemake profile. |
| `config/databases_and_containers.yaml` | Available databases, images, and runtime settings. |
| `config/datatypes.yaml`, `config/skip.yaml` | Input-format/feature declarations and default skip flags. |
| `profile/*/config.yaml` | Snakemake command-line defaults, including scheduler submission. |

Make project-specific copies before changing shared defaults. The copied core
must point to any copied resource, cluster, database, or tool-table locations.

## Precedence: core values can replace project values

The pipeline's loader changes the usual expectation that command-line values
always win:

1. Snakemake supplies project and command-line configuration.
2. Explicit project/CLI `input_dir`, `out_dir`, and `tmp_dir` values are preserved.
   For other keys, the selected core YAML replaces every colliding **top-level** value using
   `config.update`. For example, its entire `parameters` dictionary replaces a
   project dictionary of that name.
3. Secondary-tool and database settings recursively fill missing entries.
   Core-tool tables are loaded next; cluster/datatype defaults and skip flags
   follow. Existing scalar/list entries are retained by recursive default filling.
4. Resource tables assign each parameter set's resource dictionaries. Finally,
   missing entries in the selected parameter set are filled from `default`.

Put dataset paths in the project YAML, for example `out_dir: results/sample/` and
`tmp_dir: results/sample/tmp/`. Set both explicitly: overriding `out_dir` does not
automatically rewrite `tmp_dir`. Core paths are defaults when a project/CLI value
is absent. These three path keys are the exception to core precedence.
`stage_list`, `genome_prefix`, and explicit skip flags survive with the shipped
core because it does not define those keys.

A present `null` does not inherit a missing default. `input_datatypes` is a
comma-separated string, not a YAML list.

## Select stages, tools, and presets

```yaml
input_datatypes: hifi
genome_prefix: sample
stage_list: [contig]
stage_coretools:
  contig: [hifiasm, flye]
coretool_option_sets:
  hifiasm: [l0, l1_no_ul]
  flye: [hifi]
```

This requests three assembly branches: two Hifiasm presets and one Flye preset.
Each branch can introduce its own downstream analyses. Preset names must match
rows in the selected tool table. The loader converts `NA`, numeric values, and
datatype lists. Literal `default` in options handled by `parse_option` omits the
external-tool flag; it does not copy another preset. Numeric zero remains a value,
as in `-f 0`.

Secondary tools use a **single preset name**, for example:

```yaml
other_tool_option_sets:
  kmer_qc: default
  pretextmap: low_res_only
```

These selectors choose entries under `parameters.<set>.tool_options.<tool>`;
they do not add a stage by themselves.

## Override a preset in the private core

For example, disable Hifiasm's initial Bloom filter for a small genome with this
core-file fragment:

```yaml
parameters:
  micro:
    tool_options:
      hifiasm:
        l1_no_ul:
          bloom_filter_bits: 0
```

Merge this into the copied core's existing dictionary; retain its other
parameter sets and defaults. The remaining `l1_no_ul` fields are filled from
`default`. This works with the current layout, which has no `coretools/micro/`
directory. If that directory exists, table loading resets core-tool dictionaries
for `micro` before fallback. Inspect the effective configuration after changes.

## Resource profiles and parameter sets

`resource_profile: high` selects a directory; `parameter_set: micro` selects its
column and the corresponding tool-parameter set. These are separate choices.
Shipped parameter sets are `micro`, `mini`, `small`, `normal`, `big`, `large`,
`giant`, and `titan`, plus the fallback `default`.

To customize resources, copy **all three** tables from `config/resources/high/`
to a private directory such as `config/local/resources/local/`. Point the
private core's `resources_dir` to `config/local/resources/`, set
`resource_profile: local` in the project, and edit the relevant rule/column cells.
Keep table headers and rule rows: included rules can require entries even when
their tool is not selected.

`threads` controls rule threads, `memory_mb` supplies the custom `mem` resource,
and `time` supplies scheduler time strings. Values loaded from these tables
replace resource dictionaries in the core; edit the tables rather than adding a
competing `parameters.micro.threads` dictionary. See [running](running.md#resources)
for CPU-request differences and local memory limits.

## QC and dependency settings

`assembly_qc_level` selects shared assembly-QC targets per stage; `skip_*` flags
can suppress individual analyses. Set a level for every selected assembly stage.
See [stages](stages.md#quality-control-levels) for levels 0–7. Level 0 does not
remove assembly dependencies such as coverage estimation or the separate FCS
filtering path; configure their settings explicitly.

| Settings | What to check |
| --- | --- |
| `final_kmer_counter`, `final_kmer_length`, `final_kmer_datatypes` | Coverage-estimation selection; selected datatypes must exist. Initialization can add the counter/k value to active k-mer settings. |
| `data_filtering` | Eligible datatypes to filter; final-read links otherwise use raw data. |
| `ploidy`, `phasing_stage`, `assembly_kmer_length` | Modeling assumptions and read-phasing source; input presence and presets also affect assembled haplotype count. |
| `database_set` and `allowed_databases` | Requested database labels versus their actual configured paths. |
| `tool_containers`, `singularity_load_mode`, `singularity_load_env`, `singularity_load_str` | Images and runtime loading; provision the selected runtime and environment before enabling container-based tools. |
| `busco_lineage_list`, `tax_id`, `species_name`, `telomere_motif` | Sample-specific biological settings; set these for the organism being assembled. |

After parsing, initialization writes `config/config.final.yaml`,
`config/final_data.yaml`, and `requested_files.tab` under `out_dir`. Use these to
inspect effective values and targets. They are processed snapshots, not verified
replacement project YAMLs. A dry run writes these files too.

Shared QC statistics use the selected `out_dir`, including nested dataset roots.
See [summary-table limitations](outputs.md#summary-tables) for
remaining issues with metric names and thresholded totals.
