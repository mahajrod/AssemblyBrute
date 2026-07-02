ruleorder: remove_fcs_contaminants > gfa2fasta

rule fcs: # TODO: refactor to follow new scheme of the log report and new, more universal usage of the rules
    priority: 10000
    input:
        fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.unfiltered.fasta",
        db=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["path"],
        image=config["tool_containers"]["fcs_gx"],
    output:
        taxonomy=config["out_dir"] / "contig/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.contig.{haplotype}.unfiltered.{database}.taxonomy",
        summary=config["out_dir"] / "contig/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.contig.{haplotype}.unfiltered.{database}.summary"
    params:
        singularity_load_str=(config["singularity_load_str"] + "; ") if config["singularity_load_mode"] == "cluster" else "",
        tax_id=config["tax_id"]
    log:
        std=config["out_dir"] / "log/fcs.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.log",
        post=config["out_dir"] / "log/fcs.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.post.log",
        cluster_log=config["out_dir"] / "log/fcs.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.log",
        cluster_err=config["out_dir"] / "log/fcs.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fcs.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.benchmark.txt"
    conda:
        singularity_conda_env
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fcs"),
        cpus=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["threads"],
        time=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["time"],
        mem=lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["memory_mb"],
        fcs=1
    threads: lambda wildcards: config["allowed_databases"]["fcs"][wildcards.database]["threads"],

    shell: # as report(summary) might be modified manually, original version is backuped with .original extension,# || true was added as workaround to handle singularity issue with removal of rootfs after cmd
        " {params.singularity_load_str} "
        " OUT_DIR=`dirname {output.taxonomy}`; "
        " OUT_DIR=`realpath -s ${{OUT_DIR}}`; "
        " TMP_DIR=${{OUT_DIR}}'/tmp_{wildcards.database}/'; "
        " SINGULARITYENV_TMP_DIR=${{OUT_DIR}}'/singularity_{wildcards.database}/'; "
        " SINGULARITYENV_SQLITE_TMP_DIR=${{OUT_DIR}}'/singularity_sqlite_{wildcards.database}/'; "
        " mkdir -p ${{TMP_DIR}} ${{SINGULARITYENV_TMP_DIR}} ${{SINGULARITYENV_SQLITE_TMP_DIR}}; "
        " FCS_DB_DIR=`dirname {input.db}`; "
        " FCS_DB_DIR=`realpath ${{FCS_DB_DIR}}`; "
        " FASTA_DIR=`dirname {input.fasta}`; "
        " FASTA_DIR=`realpath ${{FASTA_DIR}}`; "
        " FASTA=`basename {input.fasta}`; "
        " NUM_CORES={threads}; "
        " TMPDIR=${{TMP_DIR}} SINGULARITYENV_TMPDIR=${{SINGULARITYENV_TMP_DIR}} SINGULARITYENV_SQLITE_TMPDIR=${{SINGULARITYENV_SQLITE_TMP_DIR}} "
        "     singularity exec --pid --bind ${{FCS_DB_DIR}}:/app/db/gxdb --bind ${{FASTA_DIR}}:/sample-volume/ "
        "         --bind ${{OUT_DIR}}:/output-volume/ {input.image} python3 /app/bin/run_gx --fasta /sample-volume/${{FASTA}}  "
        "         --out-dir /output-volume/ --tax-id {params.tax_id} --gx-db /app/db/gxdb/gxdb > {log.std} 2>&1 || : ; "
        " REPORT={output.taxonomy}; "
        " SUMMARY={output.summary}; "
        " cp ${{REPORT%.{wildcards.database}.taxonomy}}.{params.tax_id}.{wildcards.database}_report.txt ${{SUMMARY}}.original > {log.post} 2>&1; "
        " mv ${{REPORT%.{wildcards.database}.taxonomy}}.{params.tax_id}.{wildcards.database}_report.txt ${{SUMMARY}} >> {log.post} 2>&1; "
        " mv ${{SUMMARY%.{wildcards.database}.summary}}.{params.tax_id}.taxonomy.rpt ${{REPORT}} >> {log.post} 2>&1; "
        " rm -rf  ${{TMP_DIR}} ${{SINGULARITYENV_TMP_DIR}} ${{SINGULARITYENV_SQLITE_TMP_DIR}} >> {log.post} 2>&1; "

rule remove_fcs_contaminants: #
    priority: 5000
    input:
        fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.unfiltered.fasta",
        image=config["tool_containers"]["fcs_gx"] if not config["skip_fcs"] else [],
        fcs_report=(config["out_dir"] / ("contig/{parameters}/contamination_scan/{haplotype}/fcs/%s/{genome_prefix}.contig.{haplotype}.unfiltered.%s.summary" % (config["final_fcs_db"], config["final_fcs_db"]))) if not config["skip_fcs"] else []
    output:
        fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.fasta",
        contaminant_fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.contaminant.fasta"
    params:
        singularity_load_str=(config["singularity_load_str"] + "; ") if config["singularity_load_mode"] == "cluster" else "",
        skip="skip" if config["skip_fcs"] else "filter"
    log:
        std=config["out_dir"] / "log/remove_fcs_contaminants.contig.{parameters}.{genome_prefix}.{haplotype}.log",
        cp=config["out_dir"] / "log/remove_fcs_contaminants.contig.{parameters}.{genome_prefix}.{haplotype}.cp.log",
        rm=config["out_dir"] / "log/remove_fcs_contaminants.contig.{parameters}.{genome_prefix}.{haplotype}.rm.log",
        cluster_log=config["out_dir"] / "log/remove_fcs_contaminants.contig.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/remove_fcs_contaminants.contig.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/remove_fcs_contaminants.contig.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        singularity_conda_env
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("remove_fcs_contaminants"),
        cpus=parameters["threads"]["remove_fcs_contaminants"],
        time=parameters["time"]["remove_fcs_contaminants"],
        mem=parameters["memory_mb"]["remove_fcs_contaminants"],
        fcs=1
    threads: parameters["threads"]["remove_fcs_contaminants"],

    shell: # || true was added as workaround to handle singularity issue with removal of rootfs after cmd
        " {params.singularity_load_str} "
        " if [ '{params.skip}' = 'filter' ]; "
        " then "
        "     OUTDIR=`dirname {output.fasta}`; "
        "     OUTDIR=`realpath -s ${{OUTDIR}}`; "
        "     TMPDIR=${{OUTDIR}}'/tmp/'; "
        "     SINGULARITYENV_TMPDIR=${{OUTDIR}}'/singularity/'; "
        "     SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}'/singularity_sqlite/'; "
        "     mkdir -p ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}}; "
        "     NUM_CORES={threads}; "
        "     export FCS_DEFAULT_IMAGE={input.image}; "
        "     cat {input.fasta} | "
        "     TMPDIR=${{TMPDIR}} SINGULARITYENV_TMPDIR=${{SINGULARITYENV_TMPDIR}} SINGULARITYENV_SQLITE_TMPDIR=${{SINGULARITYENV_SQLITE_TMPDIR}} "
        "     workflow/external_tools/fcs-gx/fcs.py clean genome --action-report {input.fcs_report} "
        "     --output {output.fasta} --contam-fasta-out {output.contaminant_fasta} > {log.std} 2>&1; "
        "     rm -rf  ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}} > {log.rm} 2>&1; "
        " else "
        "     cp -f {input.fasta} {output.fasta} > {log.cp} 2>&1; "
        "     touch {output.contaminant_fasta} >> {log.cp} 2>&1; "
        " fi; "


rule fcs_adaptor: #
    priority: 2000
    input:
        fasta=config["out_dir"] / "contig/{parameters}/{genome_prefix}.contig.{haplotype}.unfiltered.fasta",
        image=config["tool_containers"]["fcs_adaptor"],
    output:
        report=config["out_dir"] / "contig/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.contig.{haplotype}.unfiltered.{database}.report",
        report_jsonl=config["out_dir"] / "contig/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.contig.{haplotype}.unfiltered.{database}.report.jsonl",
    params:
        singularity_load_str=(config["singularity_load_str"] + "; ") if config["singularity_load_mode"] == "cluster" else "",
        tax_id=config["tax_id"],
        taxonomy= lambda wildcards: " --euk " if config["allowed_databases"]["fcs_adaptor"][wildcards.database]["taxonomy"] == "eukaryota"  else " --prok "
    log:
        std=config["out_dir"] / "log/fcs_adaptor.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.log",
        post=config["out_dir"] / "log/fcs_adaptor.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.postlog",
        cluster_log=config["out_dir"] / "log/fcs_adaptor.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.cluster.log",
        cluster_err=config["out_dir"] / "log/fcs_adaptor.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.err"
    benchmark:
        config["out_dir"] / "log/fcs_adaptor.contig.{parameters}.{genome_prefix}.{haplotype}.{database}.benchmark.txt"
    conda:
        singularity_conda_env
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fcs_adaptor"),
        cpus=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["threads"],
        time=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["time"],
        mem=lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["memory_mb"],
        fcs_adaptor=1
    threads: lambda wildcards: config["allowed_databases"]["fcs_adaptor"][wildcards.database]["threads"],

    shell:
        " {params.singularity_load_str} "
        " OUTDIR=`dirname {output.report}`; "
        " OUTDIR=`realpath -s ${{OUTDIR}}`; "
        " TMPDIR=${{OUTDIR}}'/tmp/'; "
        " SINGULARITYENV_TMPDIR=${{OUTDIR}}'/singularity/'; "
        " SINGULARITYENV_SQLITE_TMPDIR=${{OUTDIR}}'/singularity_sqlite/'; "
        " mkdir -p ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}}; "
        " TMPDIR=${{TMPDIR}} SINGULARITYENV_TMPDIR=${{SINGULARITYENV_TMPDIR}} SINGULARITYENV_SQLITE_TMPDIR=${{SINGULARITYENV_SQLITE_TMPDIR}} "
        "     workflow/external_tools/fcsadaptor/run_fcsadaptor.sh --image {input.image} --fasta-input {input.fasta} "
        "         --output-dir `dirname {output.report}` {params.taxonomy} --container-engine singularity > {log.std} 2>&1; "
        " mv `dirname {output.report}`/fcs_adaptor_report.txt {output.report} > {log.post} 2>&1; "
        " mv `dirname {output.report_jsonl}`/combined.calls.jsonl {output.report_jsonl} >> {log.post} 2>&1; "
        " rm -rf ${{TMPDIR}} ${{SINGULARITYENV_TMPDIR}} ${{SINGULARITYENV_SQLITE_TMPDIR}} >> {log.post} 2>&1; "
