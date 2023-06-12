#!/bin/bash

{

#### Description: Wrapper script to phase genomic variants from ENCODE DCC hic-pipeline.
#### Usage: bash ./run-hic-phaser-encode.sh [options] <path_to_vcf> <path_to_merged_dedupped_bam_1> ... <path_to_merged_dedup_bam_N>.
#### Input: vcf file, bam file(s).
#### Output: phased vcf & phasing validation hic maps.
#### Dependencies: 3D-DNA, GNU Parallel, Java.
#### Written by: Olga Dudchenko, 01/07/2022 [updated from version dates: 09/01/2021; 08/18/2020; 10/14/2021].

echo "*****************************************************" >&1
echo "cmd log: "$0" "$* >&1
echo "*****************************************************" >&1


USAGE="
*****************************************************
Phase genomic variants after ENCODE DCC hic-pipeline.

USAGE: ./run-hic-phaser-encode.sh [options] <path_to_vcf> <path_to_merged_dedup_bam_1> ... <path_to_merged_dedup_bam_N>

DESCRIPTION:
This is a wrapper script to use 3D-DNA phaser to phase SNPs using Hi-C alignment data as generated by the ENCODE DCC hic & variant calling pipeline.

ARGUMENTS:
path_to_vcf
						Path to a Variant Call Format (vcf) file containing sequence variation data, e.g. as generated by the ENCODE DCC Hi-C variant calling pipeline.

path_to_merged_dedup_bam
						Path to bam file containing deduplicated alignments of Hi-C reads in bam format (output by Juicer2). Multiple bam files can be passed as arguments.

OPTIONS:
-h|--help
						Shows this help.

DATA FILTERING:						
-c|--chr [chr_list]
		                Phase only specific molecules (default phases all molecules that have variants listed in the vcf file). Note that -c and -C are incompatible in that if both are invoked, only the last listed option will be taken into account.

-C|--exclude-chr [chr_list]
						Remove specific molecules from the chromosome list (default: chrY). Note that -c and -C are incompatible in that if both are invoked, only the last listed option will be taken into account.

-q|--mapq [mapq]
		                Consider only Hi-C reads that align with minimal mapping quality of mapq (default is 1).

PHASER CONTROL:
-s|--stringency [stringency]
		                Specify stringency parameter for the phaser (default is 3, i.e. 3-fold enrichment in Hi-C signal to one molecule as compared to the other is necessary to phase).

-b|--background [background]
		                Specify background parameter for the phaser (default is 1, i.e. calculate enrichment on top of the noise level of 1 read)

-v|--verbose [step_for_intermediate_data_dump]
						Print intermediate phasing results every [step_for_intermediate_data_dump] steps.

OUTPUT CONTROL:
--separate-homolog-maps
						Build two separate contact maps for homologs as opposed to a single diploid contact map with interleaved homologous chromosomes. This is the preferred mode for comparing contacts across homologs using the \"Observed over Control\" view for map comparison in Juicebox. Importantly, assignment of chromosomes to homologs is arbitrary. Default: not invoked. 

WORKFLOW CONTROL:
-t|--threads [num]
        				Indicate how many threads to use. Default: half of available cores as calculated by parallel --number-of-cores.

--from-stage [pipeline_stage]
						Fast-forward to a particular stage of the pipeline. The pipeline_stage argument can be \"prep\", \"parse_vcf\", \"parse_bam\", \"visualize_input\", \"phase\", \"visualize_output\", \"update_vcf\", \"map_diploid\", \"cleanup\".

--to-stage [pipeline_stage]
						Exit after a particular stage of the pipeline. The argument can be \"prep\", \"parse_vcf\", \"parse_bam\", \"visualize_input\", \"phase\", \"visualize_output\", \"update_vcf\", \"map_diploid\", \"cleanup\".

*****************************************************
"

# path to 3D-DNA pipeline
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

# defaults:
chr=""
exclude_chr="chrY"
mapq=1
stringency=3
background=1
verbose=""
separate_homolog_maps=false

#multithreading
threads=`parallel --number-of-cores`
threads=$((threads/2))
# adjust for mem usage
tmp=`awk '/MemTotal/ {threads=int($2/1024/1024/2/6-1)}END{print threads+0}' /proc/meminfo 2>/dev/null`
tmp=$((tmp+0))
([ $tmp -gt 0 ] && [ $tmp -lt $threads ]) && threads=$tmp

#staging
first_stage="prep"
last_stage="cleanup"
declare -A stage
stage[prep]=0
stage[parse_vcf]=1
stage[parse_bam]=2
stage[visualize_input]=3
stage[phase]=4
stage[visualize_output]=5
stage[update_vcf]=6
stage[map_diploid]=7
stage[build_accessibility]=8
stage[cleanup]=9

############### HANDLE OPTIONS ###############

while :; do
	case $1 in
		-h|--help)
			echo "$USAGE" >&1
			exit 0
        ;;
### data filtration
        -c|--chr) OPTARG=$2
			echo "... -c|--chr flag was triggered, ignoring all sequences in the vcf except for $OPTARG." >&1
			chr=$OPTARG
			exclude_chr=""
        	shift
        ;;
		-C|--exclude-chr) OPTARG=$2
			echo "... -C|--exclude-chr flag was triggered, will ignore variants on $OPTARG." >&1
			exclude_chr=$OPTARG
			chr=""
        	shift
        ;;
        -q|--mapq) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -q|--mapq flag was triggered, phasing using reads with at least $OPTARG mapping quality." >&1
				mapq=$OPTARG
			else
				echo ":( Wrong syntax for mapping quality parameter. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
### phaser
        -s|--stringency) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -s|--stringency flag was triggered, phasing requiring at least $OPTARG enrichment of one haplotype vs the other." >&1
				stringency=$OPTARG
			else
				echo ":( Wrong syntax for stringency parameter. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
        -b|--background) OPTARG=$2
			re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
				echo "... -b|--background flag was triggered, phasing requiring enrichment against $OPTARG read(s) background noise." >&1
				background=$OPTARG
			else
				echo ":( Wrong syntax for background parameter. Exiting!" >&2
				exit 1
			fi
        	shift
        ;;
		-v|--verbose) OPTARG=$2
			verbose=$OPTARG
			shift
		;;
### output
		--separate-homolog-maps)
			echo "... --separate-homolog-maps flag was triggered, will build two separate contact maps (-r.hic and -a.hic) for chromosomal homologs with identical chromosomal labels." >&1
			separate_homolog_maps=true
		;;
### workflow	
        -t|--threads) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
					echo "... -t|--threads flag was triggered, will try to parallelize across $OPTARG threads." >&1
					threads=$OPTARG
			else
					echo " :( Wrong syntax for thread count parameter value. Exiting!" >&2
					exit 1
			fi        	
        	shift
        ;;	
		--from-stage) OPTARG=$2
			if [ "$OPTARG" == "prep" ] || [ "$OPTARG" == "parse_vcf" ] || [ "$OPTARG" == "parse_bam" ] || [ "$OPTARG" == "visualize_input" ] || [ "$OPTARG" == "phase" ] || [ "$OPTARG" == "visualize_output" ] || [ "$OPTARG" == "update_vcf" ] || [ "$OPTARG" == "map_diploid" ] || [ "$OPTARG" == "build_accessibility" ] || [ "$OPTARG" == "cleanup" ]; then
        		echo "... --from-stage flag was triggered. Will fast-forward to $OPTARG." >&1
        		first_stage=$OPTARG
			else
				echo " :( Whong syntax for pipeline stage. Please use prep/parse_vcf/parse_bam/visualize_input/phase/visualize_output/update_vcf/map_diploid/build_accessibility. Exiting!" >&2
				exit 1
			fi
			shift
        ;;
		--to-stage) OPTARG=$2
			if [ "$OPTARG" == "prep" ] || [ "$OPTARG" == "parse_vcf" ] || [ "$OPTARG" == "parse_bam" ] || [ "$OPTARG" == "visualize_input" ] || [ "$OPTARG" == "phase" ] || [ "$OPTARG" == "visualize_output" ] || [ "$OPTARG" == "update_vcf" ] || [ "$OPTARG" == "map_diploid" ] || [ "$OPTARG" == "build_accessibility" ] || [ "$OPTARG" == "cleanup" ]; then
				echo "... --to-stage flag was triggered. Will exit after $OPTARG." >&1
				last_stage=$OPTARG
			else
				echo " :( Whong syntax for pipeline stage. Please use prep/parse_vcf/parse_bam/visualize_input/phase/visualize_output/update_vcf/map_diploid/build_accessibility. Exiting!" >&2
				exit 1			
			fi
			shift
		;;
### utilitarian
        --) # End of all options
			shift
			break
		;;
		-?*)
			echo ":| WARNING: Unknown option. Ignoring: ${1}" >&2
		;;
		*) # Default case: If no more options then break out of the loop.
			break
	esac
	shift
done

if [[ "${stage[$first_stage]}" -gt "${stage[$last_stage]}" ]]; then
	echo >&2 ":( Please make sure that the first stage requested is in fact an earlier stage of the pipeline to the one requested as last. Exiting!"
	exit 1
fi

############### HANDLE EXTERNAL DEPENDENCIES ###############

##	Java Dependency
type java >/dev/null 2>&1 || { echo >&2 ":( Java is not available, please install/add to path Java. Exiting!"; exit 1; }

##	GNU Parallel Dependency
type parallel >/dev/null 2>&1 || { echo >&2 ":( GNU Parallel support is set to true (default) but GNU Parallel is not in the path. Please install GNU Parallel or set -p option to false. Exiting!"; exit 1; }
[ $(parallel --version | awk 'NR==1{print $3}') -ge 20150322 ] || { echo >&2 ":( Outdated version of GNU Parallel is installed. Please install/add to path v 20150322 or later. Exiting!"; exit 1; }

## Samtools Dependency
type samtools >/dev/null 2>&1 || { echo >&2 ":( Samtools are not available, please install/add to path. Exiting!"; exit 1; }
ver=`samtools --version | awk 'NR==1{print \$NF}'`
[[ $(echo "$ver < 1.13" |bc -l) -eq 1 ]] && { echo >&2 ":( Outdated version of samtools is installed. Please install/add to path v 1.13 or later. Exiting!"; exit 1; }

## TODO: kentUtils Dependency
type bedGraphToBigWig >/dev/null 2>&1 || { echo >&2 ":( bedGraphToBigWig is not available, please install/add to path, e.g. from kentUtils. Exiting!"; exit 1; }

############### HANDLE ARGUMENTS ###############

[ -z $1 ] || [ -z $2 ] && echo >&2 ":( Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE" && exit 1

[ ! -s $1 ] || [ ! -s $2 ] && echo >&2 ":( Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE" && exit 1

if [ "$#" -lt 2 ]; then
    echo >&2 "Illegal number of arguments. Please double check your input. Exiting!" && echo >&2 "$USAGE" && exit 1
fi

vcf=$1
bam=`echo "${@:2}"`
##TODO: check file extentions

## Check that reference names do not contain the \":\" character. TODO: use some more obsure separator
test=`awk -F , '$0!~/^#/{exit}($1~/##contig=<ID=/){n+=gsub(":",":",$1)}END{print n+=0}' $vcf`
[ $test -eq 0 ] || { echo ":( Sequence names contain semicolumns. This will interfere with the internal notation of the phasing pipeline. Please rename your sequences to proceed. Exiting!" | tee -a /dev/stderr && exit 1; }

############### MAIN #################

## 0. PREP BAM FILE

if [ "$first_stage" == "prep" ]; then

	echo "...Extracting unique paired alignments from bam and sorting..." >&1

	# make header for the merged file pipe
	parallel --will-cite "samtools view -H {} > {}_header.bam" ::: $bam
	header_list=`parallel --will-cite "printf %s' ' {}_header.bam" ::: $bam`
	samtools merge --no-PG -f mega_header.bam ${header_list}
	rm ${header_list}

	samtools cat -@ $((threads * 2)) -h mega_header.bam $bam | samtools view -u -d "rt:0" -d "rt:1" -d "rt:2" -d "rt:3" -d "rt:4" -d "rt:5" -@ $((threads * 2)) -F 0x400 -q $mapq - |  samtools sort -@ $threads -m 6G -o reads.sorted.bam
	[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at bam sorting. See stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }
	rm mega_header.bam

	samtools index -@ $threads reads.sorted.bam	
	[ $? -eq 0 ] || { echo ":( Failed at bam indexing. See stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }		
	# e.g. will fail with chr longer than ~500Mb. Use samtools index -c -m 14 reads.sorted.bam

	echo ":) Done extracting unique paired alignments from bam and sorting." >&1

	[ "$last_stage" == "prep" ] && { echo "Done with the requested workflow. Exiting after prepping bam!"; exit; }
	first_stage="parse_vcf"

else

	([ -f reads.sorted.bam ] && [ -f reads.sorted.bam.bai ]) || { echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr; exit 1; }

fi

## I. PARSE VCF FILE [ASSUMES THE VCF FILE IS PROPERLY SORTED]

if [ "$first_stage" == "parse_vcf" ]; then

	echo "...Parsing vcf file..." >&1
	
	awk -v chr=${chr} -v exclude_chr=${exclude_chr} -v output_prefix="in" -f ${pipeline}/phase/vcf-to-psf-and-assembly.awk ${vcf}
	mv "in.assembly" `basename ${vcf} .vcf`".in.assembly"

	assembly=`basename ${vcf} .vcf`".in.assembly"
	psf="in.psf"

	if [ -z "$chr" ]; then
		chr=$(awk '$0~/^>/{if($1!=prev){str=str"|"substr($1,2); prev=$1;}}END{print substr(str,2)}' ${psf})
	fi

	if [ ! -s $psf ] || [ ! -s `basename ${vcf} .vcf`".in.assembly" ]; then
		echo ":( Pipeline failed at parsing vcf. Check stderr for more info. Exiting!" |  tee -a /dev/stderr
		exit 1
	fi

	echo ":) Done parsing vcf file." >&1

	[ "$last_stage" == "parse_vcf" ] && { echo "Done with the requested workflow. Exiting after parsing vcf!"; exit; }
	first_stage="parse_bam"

else

	psf="in.psf"

	if [ ! -s $psf ] || [ ! -s `basename ${vcf} .vcf`".in.assembly" ]; then
		echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr
		exit 1
	fi

	if [ -z "$chr" ]; then
		chr=$(awk '$0~/^>/{if($1!=prev){str=str"|"substr($1,2); prev=$1;}}END{print substr(str,2)}' ${psf})
	fi

fi

## II. PARSE BAM FILE FOR READS OVERLAPPING SNPs.
if [ "$first_stage" == "parse_bam" ]; then
	
	echo "...Parsing bam file..." >&1

	if [[ $chr != *"|"* ]]; then
		samtools view -@ $threads reads.sorted.bam $chr | awk -f ${pipeline}/phase/extract-SNP-reads-from-sam-file.awk ${psf} - > dangling.sam

		[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at parsing bam. Check stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }
		## maybe TODO: parallelize single-chrom workflow based on interval_list
	
	else
		export SHELL=$(type -p bash)
		export psf=${psf}
		export pipeline=${pipeline}
		doit () { 
			samtools view -@ 2 reads.sorted.bam $1 | awk -f ${pipeline}/phase/extract-SNP-reads-from-sam-file.awk ${psf} -
		}
		export -f doit
		echo $chr | tr "|" "\n" | parallel -j $threads --will-cite --joblog temp.log doit > dangling.sam
		exitval=`awk 'NR>1{if($7!=0){c=1; exit}}END{print c+0}' temp.log`
		[ $exitval -eq 0 ] || { echo ":( Pipeline failed at parsing bam. Check stderr for more info. Exiting! " | tee -a /dev/stderr && exit 1; }
		rm temp.log
	fi

	awk '{c[$1]++}END{for (i in c){if(c[i]==2){print i}}}' dangling.sam | awk 'FILENAME==ARGV[1]{remember[$1]=1;next}($1 in remember){str[$1]=str[$1]" "$3}END{for(i in str){split(substr(str[i],2),a," "); print 0,a[1],1,0,0,a[2],1,1,1,"-","-",1,"-","-"}}' - dangling.sam > snp.mnd.txt
	edge_mnd="snp.mnd.txt"

	if [ ! -s ${edge_mnd} ] || [ ! -s "dangling.sam" ]; then
		echo "Pipeline failed at parsing bam. Check stderr for more info. Exiting! " | tee -a /dev/stderr
		exit 1
	fi

	echo ":) Done parsing bam file." >&1

	[ "$last_stage" == "parse_bam" ] && { echo "Done with the requested workflow. Exiting after parsing bam!"; exit; }
	first_stage="visualize_input"

else
	edge_mnd="snp.mnd.txt"

	if [ ! -s ${edge_mnd} ] || [ ! -s "dangling.sam" ]; then
		echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr
		exit 1
	fi
fi

## III. VISUALIZE INPUT PHASED BLOCKS.
if [ "$first_stage" == "visualize_input" ]; then

	echo "...Visualizing input phased blocks..." >&1

	{ bash ${pipeline}/visualize/run-assembly-visualizer.sh -c `basename ${vcf} .vcf`".in.assembly" ${edge_mnd} | sed 's/^/.../'; } 2> >(while read line; do echo "...$line" >&2; done)
 	[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Something went wrong. Check stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }

	echo ":) Done visualizing input phased blocks." >&1
	
	[ "$last_stage" == "visualize_input" ] && { echo "Done with the requested workflow. Exiting after visualizing input phased blocks!"; exit; }
	first_stage="phase"
fi

## IV. PHASE!

if [ "$first_stage" == "phase" ]; then

	echo "...Phasing..."
	if [[ $chr != *"|"* ]]; then
		{ awk -v stringency=${stringency} -v background=${background} -v outfile="out.psf" -v verbose=${verbose} -f ${pipeline}/phase/phase-intrachromosomal.awk ${psf} ${edge_mnd} | sed 's/^/.../';  } 2> >(while read line; do echo "...$line" >&2; done)
		[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at phasing. Check stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }
	else
		export SHELL=$(type -p bash)
		export psf=${psf}
		export edge_mnd=${edge_mnd}
		export stringency=${stringency}
		export background=${background}
		export verbose=${verbose}
		export pipeline=${pipeline}
		doit () { 
			cmd="echo \"Phasing chr $1.\" && awk -v chr=$1 '\$1==\">\"chr{print; id[\$NF]=1; id[-\$NF]=1}\$1~/^>/{next}(\$1 in id){print}' ${psf} > h.$1.psf && awk -v stringency=${stringency} -v background=${background} -v outfile=out.$1.psf -v verbose=${verbose} -f ${pipeline}/phase/phase-intrachromosomal.awk h.$1.psf ${edge_mnd} && rm h.$1.psf"
			eval $cmd
		}
		export -f doit
		echo $chr | tr "|" "\n" | parallel -j $threads --will-cite --joblog temp.log doit
		exitval=`awk 'NR>1{if($7!=0){c=1; exit}}END{print c+0}' temp.log`
		[ $exitval -eq 0 ] || { echo ":( Pipeline failed at phasing. See stderr for more info. Exiting! " | tee -a /dev/stderr && exit 1; }
		rm temp.log
		echo $chr | tr "|" "\n" | parallel -j $threads --will-cite -k "awk '\$0~/^>/' out.{}.psf" > out.psf
		echo $chr | tr "|" "\n" | parallel -j $threads --will-cite -k "awk '\$0!~/^>/' out.{}.psf" >> out.psf
		echo $chr | tr "|" "\n" | parallel -j $threads --will-cite rm out.{}.psf
	fi

	echo ":) Done phasing. See summary stats below:"
	echo "	------------"
	awk -f ${pipeline}/phase/print-phasing-stats-from-psf.awk out.psf
	echo "	------------"

	[ "$last_stage" == "phase" ] && { echo "Done with the requested workflow. Exiting after phasing!"; exit; }
	first_stage="visualize_output"

else
	[ -s out.psf ] || { echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr; exit 1; }
fi

## V. VISUALIZE OUTPUT PHASED BLOCKS.
if [ "$first_stage" == "visualize_output" ]; then

	echo "...Visualizing output phased blocks..." >&1

	awk -f ${pipeline}/phase/psf-to-assembly.awk "out.psf" > `basename ${vcf} .vcf`".out.assembly" 

	{ bash ${pipeline}/visualize/run-assembly-visualizer.sh -c `basename ${vcf} .vcf`".out.assembly" ${edge_mnd} | sed 's/^/.../'; } 2> >(while read line; do echo "...$line" >&2; done)

 	[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at visualizing output. Check stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }

	echo ":) Done visualizing output phased blocks." >&1
	
	[ "$last_stage" == "visualize_output" ] && { echo "Done with the requested workflow. Exiting after visualizing output phased blocks!"; exit; }
	first_stage="update_vcf"
fi

## VI. UPDATE VCF.
if [ "$first_stage" == "update_vcf" ]; then

	echo "...Updating vcf file with phasing info..." >&1

	awk -f ${pipeline}/phase/update-vcf-according-to-psf.awk out.psf ${vcf} > `basename ${vcf} .vcf`"_HiC.vcf"
	[ $? -eq 0 ] || { echo ":( Failed at updating vcf. See stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }		

	echo ":) Done updating vcf file with phasing info." >&1

	[ "$last_stage" == "update_vcf" ] && { echo "Done with the requested workflow. Exiting after updating the vcf file!"; exit; }
	first_stage="map_diploid"

fi

## VII. BUILD DIPLOID MAP FROM DANGLING READS

if [ "$first_stage" == "map_diploid" ]; then

	echo "...Building diploid contact maps from reads overlapping phased SNPs..." >&1
	
	## maybe TODO: add SNP-delimited diploid maps
	
	bash ${pipeline}/phase/assign-reads-to-homologs.sh -t ${threads} -c ${chr} out.psf dangling.sam

	samtools view -H reads.sorted.bam | grep '^@RG' | awk -F '\t' '{for(i=2;i<=NF;i++){if($i~/^ID:/){id=substr($i,4)};if($i~/^SM:/){sm=substr($i,4)}}; sub("[^a-zA-Z0-9\\.\\-]","_",sm); print id > "rg_"sm".txt"}'

while read sm
do	
		if [[ $chr != *"|"* ]]; then
			### single chromosome case

			samtools view -@ ${threads} -R "rg_"$sm".txt" -h reads.sorted.bam $chr | awk -v chr=$chr 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{if($2==chr"-r"||$2==chr"-a"){if(keep[$1]&&keep[$1]!=$2){delete keep[$1]}else{keep[$1]=$2}};next}$0~/^@SQ/{$2=$2"-r"; print; $2=substr($2,1,length($2)-2)"-a";print;next}$0~/^@/{print;next}($1 in keep)&&($7=="="||$7=="*"){$3=keep[$1];print}' reads_to_homologs.txt - | samtools sort -@ ${threads} -n -m 1G -O sam | awk '$0~/^@/{next}($1!=prev){if(n==2){sub("\t","",str); print str}; str=""; n=0}{for(i=12;i<=NF;i++){if($i~/^ip:i:/){$4=substr($i,6);break;}};str=str"\t"n"\t"$3"\t"$4"\t"n; n++; prev=$1}END{if(n==2){sub("\t","",str); print str}}' | sort -k 2,2 --parallel=${threads} -S6G > diploid.mnd.txt

			[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at building diploid contact maps. Check stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }
			## potentially TODO: parallelize single-chrom workflow based on interval_list
		else
			export SHELL=$(type -p bash)
			export psf=${psf}
			export pipeline=${pipeline}
			export sm=$sm
			doit () { 
				samtools view -@ 2  -R "rg_"$sm".txt" -h reads.sorted.bam $1 | awk -v chr=$1 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{if($2==chr"-r"||$2==chr"-a"){if(keep[$1]&&keep[$1]!=$2){delete keep[$1]}else{keep[$1]=$2}};next}$0~/^@SQ/{$2=$2"-r"; print; $2=substr($2,1,length($2)-2)"-a";print;next}$0~/^@/{print;next}($1 in keep)&&($7=="="||$7=="*"){$3=keep[$1];print}' reads_to_homologs.txt - | samtools sort -n -m 1G -O sam | awk '$0~/^@/{next}($1!=prev){if(n==2){sub("\t","",str); print str}; str=""; n=0}{for(i=12;i<=NF;i++){if($i~/^ip:i:/){$4=substr($i,6);break;}};str=str"\t"n"\t"$3"\t"$4"\t"n; n++; prev=$1}END{if(n==2){sub("\t","",str); print str}}' | sort -k 2,2 -S 6G

			}
			export -f doit
			echo $chr | tr "|" "\n" | parallel -j $threads --will-cite --joblog temp.log -k doit > diploid.mnd.txt
			exitval=`awk 'NR>1{if($7!=0){c=1; exit}}END{print c+0}' temp.log`
			[ $exitval -eq 0 ] || { echo ":( Pipeline failed at building diploid contact maps. See stderr for more info. Exiting! " | tee -a /dev/stderr && exit 1; }
			rm temp.log
		fi

		if [ "$separate_homolog_maps" == "true" ]; then
			{ awk '$2~/-r$/{gsub("-r","",$2); gsub("-r","",$6); print}' diploid.mnd.txt > tmp1.mnd.txt && bash ${pipeline}/visualize/juicebox_tools.sh pre tmp1.mnd.txt $sm"-haploid-r.hic" <(awk -v chr=$chr -F, 'BEGIN{split(chr,tmp,"|"); for(i in tmp){chrom[tmp[i]]=1}}$1!~/^#/{exit}($1!~/##contig=<ID=/){next}(substr($1,14) in chrom){split($2,a,"="); len=substr(a[2],1,length(a[2]-1)); print substr($1,14)"\t"len}' $vcf); } #&
			{ awk '$2~/-a$/{gsub("-a","",$2); gsub("-a","",$6); print}' diploid.mnd.txt > tmp2.mnd.txt && bash ${pipeline}/visualize/juicebox_tools.sh pre tmp2.mnd.txt $sm"-haploid-a.hic" <(awk -v chr=$chr -F, 'BEGIN{split(chr,tmp,"|"); for(i in tmp){chrom[tmp[i]]=1}}$1!~/^#/{exit}($1!~/##contig=<ID=/){next}(substr($1,14) in chrom){split($2,a,"="); len=substr(a[2],1,length(a[2]-1)); print substr($1,14)"\t"len}' $vcf); } #&
			#wait
			rm tmp1.mnd.txt tmp2.mnd.txt
			## TODO: check if successful
		else
			bash ${pipeline}/visualize/juicebox_tools.sh pre diploid.mnd.txt $sm"-diploid.hic" <(awk -v chr=$chr -F, 'BEGIN{split(chr,tmp,"|"); for(i in tmp){chrom[tmp[i]]=1}}$1!~/^#/{exit}($1!~/##contig=<ID=/){next}(substr($1,14) in chrom){split($2,a,"="); len=substr(a[2],1,length(a[2]-1)); print substr($1,14)"-r\t"len; print substr($1,14)"-a\t"len}' $vcf)
			## TODO: check if successful
		fi

	rm "rg_"$sm".txt" diploid.mnd.txt

	done < <(samtools view -H reads.sorted.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | awk '{gsub("[^a-zA-Z0-9\\.\\-]", "_")}1' | uniq)

	echo ":) Done building diploid contact maps from reads overlapping phased SNPs." >&1

	[ "$last_stage" == "map_diploid" ] && { echo "Done with the requested workflow. Exiting after building diploid hic maps!"; exit; }
	first_stage="build_accessibility"

fi

## VIII. BUILD ACCESSIBILITY TRACKS

if [ "$first_stage" == "build_accessibility" ]; then

	echo "...Building diploid accessibility tracks from reads overlapping phased SNPs..." >&1

	if [ ! -s reads.sorted.bam ] || [ ! -s reads_to_homologs.txt ]; then
		echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr
		exit 1
	fi
	
	cmd="samtools view -H reads.sorted.bam | grep '^@RG' | awk -F '\t' '{for(i=2;i<=NF;i++){if(\$i~/^ID:/){id=substr(\$i,4)};if(\$i~/^SM:/){sm=substr(\$i,4)};if(\$i~/^PL:/){print substr(\$i,4)}}; sub(\"[^a-zA-Z0-9\\\.\\\-]\",\"_\",sm); print id > \"rg_\"sm\".txt\"}' | uniq | xargs" && pl=`eval $cmd`

	# check that platforms arent' mixed. TODO: group readgroups by sm and platform
	([ "$pl" == "ILLUMINA" ] || [ "$pl" == "LS454" ]) || { echo ":( Data from different platforms seems to be mixed. Can't handle this case. Exiting!" | tee -a /dev/stderr && exit 1; }

	[ "$pl" == "ILLUMINA" ] && junction_rt_string="-d rt:2 -d rt:3 -d rt:4 -d rt:5" || junction_rt_string="-d rt:0 -d rt:1"

	while read sm
	do	
		if [[ $chr != *"|"* ]]; then
			### single chromosome case

			samtools view -@ ${threads} -R "rg_"$sm".txt" ${junction_rt_string} -h reads.sorted.bam $chr | awk -v chr=$chr 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{if($2==chr"-r"||$2==chr"-a"){if(keep[$1]&&keep[$1]!=$2){delete keep[$1]}else{keep[$1]=$2}};next}$0~/^@/{next}($1 in keep){$3=keep[$1]; for (i=12; i<=NF; i++) {if ($i~/^ip/) {split($i, ip, ":"); locus[$3" "ip[3]]++; break}}}END{for (i in locus) {split(i, a, " "); print a[1], a[2]-1, a[2], locus[i]}}' reads_to_homologs.txt - | sort -k1,1 -k2,2n --parallel=${threads} -S 6G > tmp.bedgraph

			[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at building diploid accessibility tracks. Check stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }
		else
			export SHELL=$(type -p bash)
			export psf=${psf}
			export pipeline=${pipeline}
			export sm=$sm
			export junction_rt_string=${junction_rt_string}
			doit () { 
				samtools view -@ 2  -R "rg_"$sm".txt" ${junction_rt_string} -h reads.sorted.bam $1 | awk -v chr=$1 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{if($2==chr"-r"||$2==chr"-a"){if(keep[$1]&&keep[$1]!=$2){delete keep[$1]}else{keep[$1]=$2}};next}$0~/^@/{next}($1 in keep){$3=keep[$1]; for (i=12; i<=NF; i++) {if ($i~/^ip/) {split($i, ip, ":"); locus[$3" "ip[3]]++; break}}}END{for (i in locus) {split(i, a, " "); print a[1], a[2]-1, a[2], locus[i]}}' reads_to_homologs.txt - | sort -k1,1 -k2,2n -S 6G
			}
			export -f doit
			echo $chr | tr "|" "\n" | parallel -j $threads --will-cite --joblog temp.log -k doit > tmp.bedgraph
			exitval=`awk 'NR>1{if($7!=0){c=1; exit}}END{print c+0}' temp.log`
			[ $exitval -eq 0 ] || { echo ":( Pipeline failed at building diploid contact maps. See stderr for more info. Exiting! " | tee -a /dev/stderr && exit 1; }
			rm temp.log
		fi

		bedGraphToBigWig tmp.bedgraph <(awk -v chr=$chr -F, 'BEGIN{split(chr,tmp,"|"); for(i in tmp){chrom[tmp[i]]=1}}$1!~/^#/{exit}($1!~/##contig=<ID=/){next}(substr($1,14) in chrom){split($2,a,"="); len=substr(a[2],1,length(a[2]-1)); print substr($1,14)"-r\t"len; print substr($1,14)"-a\t"len}' $vcf) $sm"-diploid.bw"

		rm "rg_"$sm".txt" tmp.bedgraph

	done < <(samtools view -H reads.sorted.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | awk '{gsub("[^a-zA-Z0-9\\.\\-]", "_")}1' | uniq)

	echo ":) Done building diploid accessibility tracks from reads overlapping phased SNPs." >&1

	[ "$last_stage" == "build_accessibility" ] && { echo "Done with the requested workflow. Exiting after building diploid accessibility tracks!"; exit; }
	first_stage="cleanup"

fi

## IX. CLEANUP
	echo "...Starting cleanup..." >&1
	#rm reads.sorted.bam reads.sorted.bam.bai
	#rm dangling.sam
	#rm snp.mnd.txt reads_to_homologs.txt
	#rm *_track.txt
	echo ":) Done with cleanup. This is the last stage of the pipeline. Exiting!"
	exit
}