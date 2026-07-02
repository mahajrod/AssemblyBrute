ruleorder:  bam2bed_for_hic_map > bam2bed
rule bam2bed_for_hic_map: # TODO: move code to bam2bed, do better management for resources
    input:
        bam="{bam_dir}/{bam_prefix}.rmdup.bam",
        log="{bam_dir}/log/",
    output:
        bed="{bam_dir}/{bam_prefix}.rmdup.bed"
    log:
        samtools="{bam_dir}/log/bam2bed_for_hic_map.{bam_prefix}.samtools.log",
        bam2bed="{bam_dir}/log/bam2bed_for_hic_map.{bam_prefix}.bam2bed.log",
        sort="{bam_dir}/log/bam2bed_for_hic_map.{bam_prefix}.sort.log",
        cluster_log="{bam_dir}/logbam2bed_for_hic_map.{bam_prefix}.cluster.log",
        cluster_err="{bam_dir}/log/bam2bed_for_hic_map.{bam_prefix}.cluster.err"
    benchmark:
        "{bam_dir}/log/bam2bed_for_hic_map.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bam2bed_for_hic_map"),
        cpus=parameters["threads"]["bam2bed_for_hic_map"],
        time=parameters["time"]["bam2bed_for_hic_map"],
        mem=parameters["memory_mb"]["bam2bed_for_hic_map"]
    threads: parameters["threads"]["bam2bed_for_hic_map"]
    shell:
        " TMP_DIR={output.bed}.tmp; "
        " mkdir -p ${{TMP_DIR}}; "
        " samtools view -@ 4 -u -F0x400 {input.bam} 2>{log.samtools} | bamToBed 2>{log.bam2bed} | "
        "     sort -k4 --parallel=20 -S{resources.mem}M -T ${{TMP_DIR}} > {output.bed} 2>{log.sort}; "
        " rm -r ${{TMP_DIR}}; "

rule bed2pairs_for_hic_map:
    input:
        bed="{bam_dir}/{bam_prefix}.rmdup.bed",
        log="{bam_dir}/log/",
    output:
        pairs="{bam_dir}/{bam_prefix}.rmdup.pairs",
    log:
        paste="{bam_dir}/log/bed2pairs.{bam_prefix}.paste.log",
        awk1="{bam_dir}/log/bed2pairs.{bam_prefix}.awk1.log",
        awk2="{bam_dir}/log/bed2pairs.{bam_prefix}.awk2.log",
        tr="{bam_dir}/log/bed2pairs.{bam_prefix}.tr.log",
        sort="{bam_dir}/log/bed2pairs.{bam_prefix}.sort.log",
        cluster_log="{bam_dir}/log/bed2pairs.{bam_prefix}.cluster.log",
        cluster_err="{bam_dir}/log/bed2pairs.{bam_prefix}.cluster.err"
    benchmark:
        "{bam_dir}/log/bed2pairs.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bed2pairs"),
        cpus=parameters["threads"]["bed2pairs"],
        time=parameters["time"]["bed2pairs"],
        mem=parameters["memory_mb"]["bed2pairs"]
    threads: parameters["threads"]["bed2pairs"]
    shell:
        " TMP_DIR={output.pairs}.tmp; "
        " mkdir -p ${{TMP_DIR}}; "
        " paste -d '\t' - - < {input.bed} 2>{log.paste} | "
        " awk 'BEGIN {{FS=\"\t\"; OFS=\"\t\"}} "
        "     {{if ($1 > $7) {{print substr($4,1,length($4)-2),$12,$7,$8,\"16\",$6,$1,$2,\"8\",$11,$5}} "
        "     else {{print substr($4,1,length($4)-2),$6,$1,$2,\"8\",$12,$7,$8,\"16\",$5,$11}} }}' 2>{log.awk1} | "
        " tr '\-+' '01' 2>{log.tr} | sort -k3,3d -k7,7d --parallel=20 -S{resources.mem}M -T ${{TMP_DIR}} 2>{log.sort} | awk 'NF==11' > {output.pairs} 2>{log.awk2};"
        " rm -r ${{TMP_DIR}}; "

rule create_higlass_track_from_bed: #
    input:
        bam="{bam_dir}/{bam_prefix}.rmdup.bam",
        csi="{bam_dir}/{bam_prefix}.rmdup.bam.csi",
        pairs="{bam_dir}/{bam_prefix}.rmdup.pairs",
        log="{bam_dir}/log/",
    output:
        genome_higlass="{bam_dir}/{bam_prefix}.higlass.genome",
        higlass_cool="{bam_dir}/{bam_prefix}.higlass.cool",
        higlass_mcool="{bam_dir}/{bam_prefix}.higlass.mcool",
    log:
        idxstat="{bam_dir}/log/create_higlass_track.{bam_prefix}.idxstat.log",
        cut="{bam_dir}/log/create_higlass_track.{bam_prefix}.cut.log",
        awk="{bam_dir}/log/create_higlass_track.{bam_prefix}.awk.log",
        sed1="{bam_dir}/log/create_higlass_track.{bam_prefix}.sed1.log",
        sort="{bam_dir}/log/create_higlass_track.{bam_prefix}.sort.log",
        cload="{bam_dir}/log/create_higlass_track.{bam_prefix}.cload.log",
        zoomify="{bam_dir}/log/create_higlass_track.{bam_prefix}.zoomify.log",
        cluster_log="{bam_dir}/log/create_higlass_track.{bam_prefix}.cluster.log",
        cluster_err="{bam_dir}/log/create_higlass_track.{bam_prefix}.cluster.err"
    benchmark:
        "{bam_dir}/log/create_higlass_track.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_higlass_track_from_bed"),
        cpus=parameters["threads"]["create_higlass_track"] ,
        time=parameters["time"]["create_higlass_track"],
        mem=parameters["memory_mb"]["create_higlass_track"]
    threads: parameters["threads"]["create_higlass_track"]
    shell:
        " samtools idxstat {input.bam} 2>{log.idxstat} | cut -f1,2 2>{log.cut} | "
        "     awk 'NR>1 {{print PREV_LINE}}; {{PREV_LINE=$0}}; END {{if (PREV_LINE !~ /^*/) {{print PREV_LINE}} }}' 2>{log.awk} | "
        "     sed 's/-/_/g' 2>{log.sed1} | sort -k2,2 -nr > {output.genome_higlass} 2>{log.sort}; "
        " cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 {output.genome_higlass}:1000 {input.pairs} {output.higlass_cool} 2>{log.cload}; "
        " cooler zoomify --resolutions 5000,10000,25000,50000,100000,150000,200000,300000,400000,500000,1000000,2500000 "
        "     -o {output.higlass_mcool} {output.higlass_cool} 2>{log.zoomify}; "
