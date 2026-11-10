# Troubleshooting

## Initialization and input problems

| Symptom | Check and action |
| --- | --- |
| Outputs appear under an unexpected root | Explicit project/CLI `input_dir`, `out_dir`, and `tmp_dir` override core defaults; other colliding core keys still take precedence. Inspect `config/config.final.yaml`. |
| Selected datatype has no input files | Discovery uses `<input_dir>/<datatype>/<format>/` without recursion. Check the comma-separated `input_datatypes` and [supported formats](inputs.md#read-directories-and-formats). |
| Mixed formats/extensions rejected | Use one recognized format and extension per datatype. |
| Paired-read count or suffix error | Files are paired after sorting; use an even count and mate names differing by one character. See [paired-read names](inputs.md#paired-read-names). |
| Draft count/haplotype suffix error | Match the uncompressed `.fasta` count to `ploidy`; multi-haplotype names must end in `hap1.fasta`, `hap2.fasta`, etc. |
| `KeyError: 'hifi'` with draft-only input | Included MitoHiFi rules access HiFi during parsing even when mtDNA is disabled. Retain HiFi inputs for draft QC. |
| Missing reference directory or required file | Use `reference/<name>/` with one `.fasta`, `.whitelist`, and `.orderlist`; the error identifies the missing type and directory. Supply `.syn` explicitly to avoid automatic creation. |
| `KeyError` for a preset or tool | Check tool-table row names and project mappings against [tool support](stages.md#tool-support). |

## Empty or surprising job graphs

| Symptom | Explanation and action |
| --- | --- |
| Selected stage requests no work | Check `requested_files.tab`; some stages have no targets without a downstream consumer. |
| `filter_reads` requests no files | Check `data_filtering` and the datatype's filtering capability. |
| No alignment BAM after adding `hic_alignment` | Select `hic_scaffolding` or an explicit BAM target; `hic_alignment` adds no BAM target to `rule all`. |
| No read-phasing outputs | Check that `phasing_stage` has multiple assembled haplotypes and eligible reads. |
| Unexpectedly few assembly branches | Check stage ordering, inferred assembly ploidy, and compatible aligner/scaffolder pairs. See [stages](stages.md). |
| Additional ploidy models do not appear | Use `ploidy_test_list` and `kmer_test_list`; the older `kmer_ploidy_test_list` key is ignored. |
| Extra dependencies despite QC level 0 | Hifiasm coverage estimation and FCS filtering are outside the shared-QC level gate. |

## Environments, resources, and reports

| Symptom | Check and action |
| --- | --- |
| Conda YAML cannot be found | `conda.<tool>.yaml` resolves relative to the main `Snakefile`. See [environment settings](installation.md#create-tool-environments). |
| Small-genome Hifiasm startup uses excessive RAM | Set `bloom_filter_bits: 0` in a private preset; see [preset overrides](configuration.md#override-a-preset-in-the-private-core). |
| Tool threads differ from scheduler CPU requests | Rule `threads` and `resources.cpus` are separate; inspect both and the [resource tables](running.md#resources). |
| Statistics columns have `.1` suffixes | The gathering script emits duplicate labels; select exact headers and compare with the source QUAST report. |
| Statistics argument error with an empty Merqury selection | The command can emit bare `-m`, even at QC level 1. Retain a valid datatype selection or fix command construction. |
| No BUSCO/Merqury/Pretext report | Check shared QC level, skip flags, and required lineages, reads, or tracks. |
| Polyploid Merqury omits extra haplotypes | The integration passes only `hap0` or a two-haplotype pair. Evaluate additional haplotypes separately. |

For incomplete integrations, see [tool support](stages.md#tool-support).
