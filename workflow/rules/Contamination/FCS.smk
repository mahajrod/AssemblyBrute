ruleorder: remove_fcs_contaminants > gfa2fasta

rule fcs: #
    priority: 10000
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.fasta",
        db=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["path"],
        image=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["image_path"],
    output:
        taxonomy=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.{database}.taxonomy",
        summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.{database}.summary"
        #report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.taxonomy.txt",
        #summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.fcs_gx_report.txt"
    params:
        tax_id=config["tax_id"]
    log:
        std=output_dict["log"]  / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.log",
        post=output_dict["log"]  / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.post.log",
        cluster_log=output_dict["cluster_log"] / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "fcs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        cpus=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["threads"],
        time=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["time"],
        mem=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["memory_mb"],
        fcs=1
    threads: lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["threads"],

    shell: # as report(summary) might be modified manually, original version is backuped with .original extension,# || true was added as workaround to handle singularity issue with removal of rootfs after cmd
        " OUTDIR=`dirname {output.taxonomy}`; "
        " OUTDIR=`realpath -s ${{OUTDIR}}`; "
        " TMPDIR=${{OUTDIR}}'/tmp_{wildcards.database}/'; "
        " SINGULARITYENV_TMPDIR=${{OUTDIR}}'/singularity_{wildcards.database}/'; "
        " SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}'/singularity_sqlite_{wildcards.database}/'; "
        " mkdir -p ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}}; "
        " NUM_CORES={threads}; "
        " export FCS_DEFAULT_IMAGE={input.image}; "
        " TMPDIR=${{TMPDIR}} SINGULARITYENV_TMPDIR=${{SINGULARITYENV_TMPDIR}} SINGULARITYENV_SQLITE_TMPDIR=${{SINGULARITYENV_SQLITE_TMPDIR}} "
        " workflow/external_tools/fcs-gx/fcs.py  screen genome --fasta {input.fasta} --out-dir `dirname {output.taxonomy}` --tax-id {params.tax_id} --gx-db {input.db} > {log.std} 2>&1 || : ; "
        " REPORT={output.taxonomy}; "
        " SUMMARY={output.summary}; "
        " cp ${{REPORT%.{wildcards.database}.taxonomy}}.{params.tax_id}.{wildcards.database}_report.txt ${{SUMMARY}}.original > {log.post} 2>&1; "
        " mv ${{REPORT%.{wildcards.database}.taxonomy}}.{params.tax_id}.{wildcards.database}_report.txt ${{SUMMARY}} >> {log.post} 2>&1; "
        " mv ${{SUMMARY%.{wildcards.database}.summary}}.{params.tax_id}.taxonomy.rpt ${{REPORT}} >> {log.post} 2>&1; "
        " rm -rf  ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}} >> {log.post} 2>&1; "

rule remove_fcs_contaminants: #
    priority: 5000
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.fasta",
        image=lambda wildcards: config["allowed_databases"]["fcs"][config["final_fcs_db"]]["image_path"],
        fcs_report=(out_dir_path / ("{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/%s/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.%s.summary" % (config["final_fcs_db"], config["final_fcs_db"]))) if not config["skip_fcs"] else []
    output:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        contaminant_fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.contaminant.fasta"
    params:
        skip="skip" if config["skip_fcs"] else "filter"
    log:
        std=output_dict["log"]  / "remove_fcs_contaminants.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cp=output_dict["log"]  / "remove_fcs_contaminants.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cp.log",
        rm=output_dict["log"]  / "remove_fcs_contaminants.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.rm.log",
        cluster_log=output_dict["cluster_log"] / "remove_fcs_contaminants.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "remove_fcs_contaminants.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "remove_fcs_contaminants.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        cpus=parameters["threads"]["remove_fcs_contaminants"],
        time=parameters["time"]["remove_fcs_contaminants"],
        mem=parameters["memory_mb"]["remove_fcs_contaminants"],
        fcs=1
    threads: parameters["threads"]["remove_fcs_contaminants"],

    shell: # || true was added as workaround to handle singularity issue with removal of rootfs after cmd
        " if [ '{params.skip}' = 'filter' ]; "
        " then "
        "       OUTDIR=`dirname {output.fasta}`; "
        "       OUTDIR=`realpath -s ${{OUTDIR}}`; "
        "       TMPDIR=${{OUTDIR}}'/tmp/'; "
        "       SINGULARITYENV_TMPDIR=${{OUTDIR}}'/singularity/'; "
        "       SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}'/singularity_sqlite/'; "
        "       mkdir -p ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}}; "
        "       NUM_CORES={threads}; "
        "       export FCS_DEFAULT_IMAGE={input.image}; "
        "       cat {input.fasta} | "
        "       TMPDIR=${{TMPDIR}} SINGULARITYENV_TMPDIR=${{SINGULARITYENV_TMPDIR}} SINGULARITYENV_SQLITE_TMPDIR=${{SINGULARITYENV_SQLITE_TMPDIR}} "
        "       workflow/external_tools/fcs-gx/fcs.py clean genome --action-report {input.fcs_report} "
        "       --output {output.fasta} --contam-fasta-out {output.contaminant_fasta} > {log.std} 2>&1; "
        "       rm -rf  ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}} > {log.rm} 2>&1; "
        " else "
        "       cp -f {input.fasta} {output.fasta} > {log.cp} 2>&1; "
        "       touch {output.contaminant_fasta} >> {log.cp} 2>&1; "
        " fi; "


rule fcs_adaptor: #
    priority: 2000
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.fasta",
        #db=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["path"],
        image=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["image_path"],
    output:
        report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.report",
        report_jsonl=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.report.jsonl",
        #summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{database}.summary"
        #report=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.taxonomy.txt",
        #summary=out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.{tax_id}.fcs_gx_report.txt"
    params:
        tax_id=config["tax_id"],
        taxonomy= lambda wildcards: " --euk " if config["allowed_databases"]["fcs_adaptor"][wildcards.database]["taxonomy"] == "eukaryota"  else " --prok "
    log:
        std=output_dict["log"]  / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.log",
        post=output_dict["log"]  / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.postlog",
        cluster_log=output_dict["cluster_log"] / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.err"
    benchmark:
        output_dict["benchmark"]  / "fcs_adaptor.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{database}.benchmark.txt"
    conda:
        config["conda"]["singularity"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["singularity"]["yaml"])
    resources:
        cpus=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["threads"],
        time=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["time"],
        mem=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["memory_mb"],
        fcs_adaptor=1
    threads: lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["threads"],

    shell:
        " OUTDIR=`dirname {output.report}`; "
        " OUTDIR=`realpath -s ${{OUTDIR}}`; "
        " TMPDIR=${{OUTDIR}}'/tmp/'; "
        " SINGULARITYENV_TMPDIR=${{OUTDIR}}'/singularity/'; "
        " SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}'/singularity_sqlite/'; "
        " mkdir -p ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}}; "
        " TMPDIR=${{TMPDIR}} SINGULARITYENV_TMPDIR=${{SINGULARITYENV_TMPDIR}} SINGULARITYENV_SQLITE_TMPDIR=${{SINGULARITYENV_SQLITE_TMPDIR}} "
        " workflow/external_tools/fcsadaptor/run_fcsadaptor.sh --image {input.image} --fasta-input {input.fasta} --output-dir `dirname {output.report}` {params.taxonomy} --container-engine singularity > {log.std} 2>&1; "
        " mv `dirname {output.report}`/fcs_adaptor_report.txt {output.report} > {log.post} 2>&1; "
        " mv `dirname {output.report_jsonl}`/combined.calls.jsonl {output.report_jsonl} >> {log.post} 2>&1; "
        " rm -rf ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}} >> {log.post} 2>&1; "
