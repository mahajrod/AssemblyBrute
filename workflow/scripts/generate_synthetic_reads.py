#!/usr/bin/env python3
"""Generate self-contained, reproducible read fixtures; Hi-C contacts are a toy model."""
import argparse
from bisect import bisect_left
from collections import Counter
from contextlib import ExitStack
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import tempfile


def format_fasta(records):
    return "".join(f">{name}\n" + "".join(sequence[i:i + 60] + "\n" for i in range(0, len(sequence), 60))
                   for name, sequence in records.items())


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chromosome-lengths", type=int, nargs=2, default=[750000, 250000],
                        help="Lengths of chr1 and chr2 (default: 750000 250000)")
    parser.add_argument("--gc-content", type=float, default=0.4,
                        help="GC fraction of the randomly generated chromosomes (default: 0.4)")
    parser.add_argument("--heterozygosity", type=float, default=0.01,
                        help="SNP fraction in a second haplotype (default: 0.01); 0 selects one haplotype")
    args = parser.parse_args()
    if not 0 <= args.heterozygosity <= 1:
        parser.error("--heterozygosity must be between 0 and 1")
    if any(length < 1000 for length in args.chromosome_lengths):
        parser.error("Each chromosome must contain at least 1000 bases")
    if not 0 < args.gc_content < 1:
        parser.error("--gc-content must be between 0 and 1, excluding the endpoints")
    reference_rng = random.Random(args.seed + 1)
    weights = [(1 - args.gc_content) / 2, args.gc_content / 2,
               args.gc_content / 2, (1 - args.gc_content) / 2]
    source_records = {f"chr{index}": "".join(reference_rng.choices("ACGT", weights=weights, k=length))
                      for index, length in enumerate(args.chromosome_lengths, 1)}
    reference_bytes = format_fasta(source_records).encode("ascii")
    records = source_records
    origins = {name: (name, 1) for name in records}
    variants = []
    ploidy = 2 if args.heterozygosity else 1
    if ploidy == 2:
        records, origins = {}, {}
        variant_rng = random.Random(args.seed)
        for name, sequence in source_records.items():
            mutated = list(sequence)
            positions = sorted(variant_rng.sample(range(len(sequence)), round(len(sequence) * args.heterozygosity)))
            for position in positions:
                alternate = variant_rng.choice([base for base in "ACGT" if base != sequence[position]])
                mutated[position] = alternate
                variants.append((name, position, sequence[position], alternate))
            for haplotype, bases in [(1, sequence), (2, "".join(mutated))]:
                label = f"{name}.hap{haplotype}"
                records[label] = bases
                origins[label] = (name, haplotype)
        if not variants:
            parser.error("--heterozygosity produces no SNPs on this reference")
    reference_length = sum(len(sequence) for sequence in source_records.values())
    rng = random.Random(args.seed)
    names = list(records)
    lengths = [len(records[name]) for name in names]
    anchors = {}
    for name, sequence in records.items():
        anchors[name] = []
        for cut in range(len(sequence) - 3):
            if sequence[cut:cut + 4] == "GATC":
                if cut + 150 <= len(sequence):
                    anchors[name].append((cut, "+"))
                if cut + 4 >= 150:
                    anchors[name].append((cut + 4 - 150, "-"))
        if not anchors[name]:
            raise ValueError(f"No usable GATC anchors for Hi-C in {name}")
        anchors[name].sort()
    anchor_positions = {name: [start for start, _ in positions] for name, positions in anchors.items()}
    hic_pairs = max(1000, math.ceil(reference_length * 36 / 300))

    outputs = {
        "hifi": "hifi/fastq/synthetic.fastq.gz",
        "wgs1": "illumina/fastq/synthetic_1.fastq.gz",
        "wgs2": "illumina/fastq/synthetic_2.fastq.gz",
        "hic1": "hic/fastq/synthetic_1.fastq.gz",
        "hic2": "hic/fastq/synthetic_2.fastq.gz",
        "truth": "truth.tsv.gz",
    }
    draft_paths = [f"draft/fasta/synthetic.hap{haplotype}.fasta" for haplotype in [1, 2]] if ploidy == 2 else ["draft/fasta/synthetic.hap0.fasta"]
    extra = ["synthetic.fasta", "manifest.json", *draft_paths]
    if ploidy == 2:
        extra.append("variants.tsv")
    directories = {".", "hifi", "hifi/fastq", "illumina", "illumina/fastq",
                   "hic", "hic/fastq", "draft", "draft/fasta"}
    markers = [str(Path(directory) / "DESCRIPTION") for directory in sorted(directories)
               if not (args.output / directory / "DESCRIPTION").exists()
               and not (args.output / directory / "DESCRIPTION").is_symlink()]
    for relative in [*outputs.values(), *extra, *markers]:
        destination = args.output / relative
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(f"Refusing to overwrite {destination}")

    args.output.mkdir(parents=True, exist_ok=True)
    counts, bases, contacts = Counter(), Counter(), Counter()
    haplotype_bases = {haplotype: Counter() for haplotype in range(1, ploidy + 1)}
    with tempfile.TemporaryDirectory(prefix=".generating-", dir=args.output) as temporary:
        stage = Path(temporary)
        with ExitStack() as stack:
            handles = {}
            for key, relative in outputs.items():
                path = stage / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                raw = stack.enter_context(path.open("xb"))
                compressed = stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0))
                handles[key] = stack.enter_context(io.TextIOWrapper(compressed, encoding="ascii", newline="\n"))
            handles["truth"].write("read_id\tdatatype\tmate\tcontig\tstart\tend\tstrand\tsubstitutions\tcontact\n")

            def emit(key, identifier, datatype, mate, contig, start, length, strand, contact="."):
                sequence = records[contig][start:start + length]
                if strand == "-":
                    sequence = reverse_complement(sequence)
                substitutions = 0
                mutated = []
                for base in sequence:
                    if rng.random() < 0.001:
                        base = rng.choice([b for b in "ACGT" if b != base])
                        substitutions += 1
                    mutated.append(base)
                handles[key].write(f"@{identifier}\n{''.join(mutated)}\n+\n{'?' * length}\n")
                handles["truth"].write(f"{identifier}\t{datatype}\t{mate}\t{contig}\t{start}\t{start + length}\t{strand}\t{substitutions}\t{contact}\n")
                counts[datatype] += 1
                bases[datatype] += length
                haplotype_bases[origins[contig][1]][datatype] += length

            for contig, reference in records.items():
                generated, count = 0, 0
                while generated < 40 / ploidy * len(reference):
                    length = rng.randint(1000, min(3000, len(reference)))
                    start = 0 if count == 0 else len(reference) - length if count == 1 else rng.randrange(len(reference) - length + 1)
                    emit("hifi", f"synthetic/{counts['hifi'] + 1}/ccs", "hifi", 0,
                         contig, start, length, rng.choice(["+", "-"]))
                    generated += length
                    count += 1
                for _ in range(math.ceil(50 / ploidy * len(reference) / 300)):
                    insert = max(300, min(len(reference), round(rng.gauss(350, 30))))
                    start = rng.randrange(len(reference) - insert + 1)
                    identifier = f"synthetic_wgs_{counts['illumina'] // 2 + 1}"
                    ends = [(start, "+"), (start + insert - 150, "-")]
                    if rng.choice([True, False]):
                        ends.reverse()
                    for mate, (position, strand) in enumerate(ends, 1):
                        emit(f"wgs{mate}", f"{identifier}/{mate}", "illumina", mate,
                             contig, position, 150, strand)

            for pair in range(hic_pairs):
                first = rng.choices(names, weights=lengths)[0]
                others = [name for name in names if name != first and origins[name][1] == origins[first][1]]
                second = rng.choice(others) if others and rng.random() < 0.1 else first
                contact = "cis" if first == second else "trans"
                contacts[contact] += 1
                start, strand = rng.choice(anchors[first])
                candidates = anchors[second]
                if contact == "cis":
                    direction = rng.choice([-1, 1])
                    available = start if direction < 0 else len(records[first]) - 150 - start
                    distance = 150 * math.expm1(rng.random() * math.log1p(available / 150))
                    target = start + direction * distance
                    index = bisect_left(anchor_positions[second], target)
                    index = min(index, len(candidates) - 1)
                    if index and abs(candidates[index - 1][0] - target) < abs(candidates[index][0] - target):
                        index -= 1
                    other_start, other_strand = candidates[index]
                else:
                    other_start, other_strand = rng.choice(candidates)
                emit("hic1", f"synthetic_hic_{pair + 1}/1", "hic", 1, first, start, 150, strand, contact)
                emit("hic2", f"synthetic_hic_{pair + 1}/2", "hic", 2, second, other_start, 150, other_strand, contact)

        for haplotype, relative in enumerate(draft_paths, 1):
            draft = stage / relative
            draft.parent.mkdir(parents=True, exist_ok=True)
            if ploidy == 1:
                draft.write_bytes(reference_bytes)
            else:
                draft.write_text(format_fasta({name: sequence for name, sequence in records.items()
                                               if origins[name][1] == haplotype}))
        if ploidy == 2:
            (stage / "variants.tsv").write_text("contig\tposition\tref\talt\n" + "".join(
                f"{name}\t{position}\t{ref}\t{alt}\n" for name, position, ref, alt in variants))
        (stage / "synthetic.fasta").write_bytes(reference_bytes)
        manifest = {
            "seed": args.seed,
            "reference_sha256": hashlib.sha256(reference_bytes).hexdigest(),
            "reference_lengths": {name: len(sequence) for name, sequence in source_records.items()},
            "models": {"hifi": "40x, 1000-3000 bp capped by contig length" if ploidy == 1 else "40x total, 20x per haplotype; 1000-3000 bp capped by contig length",
                       "illumina": "50x, paired 150 bp, FR, insert 350 +/- 30 bp (minimum 300)" if ploidy == 1 else "50x total, 25x per haplotype; paired 150 bp, FR, insert 350 +/- 30 bp (minimum 300)",
                       "hic": f"{hic_pairs} pairs, paired 150 bp, GATC anchors; toy cis distance bias, 10% trans",
                       "errors": "0.001 substitution probability, constant Q30; no indels/adapters"},
            "reads": dict(counts), "bases": dict(bases), "hic_contacts": dict(contacts),
            "coverage": {key: value / reference_length for key, value in bases.items()},
            "truth_coordinates": "0-based, half-open; substitutions counted after strand orientation",
            "sha256": {str(p.relative_to(stage)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(stage.rglob("*")) if p.is_file()},
        }
        manifest["reference_generation"] = {"seed": args.seed + 1, "gc_content": args.gc_content,
                                             "method": "Independent A/C/G/T bases; A=T and C=G"}
        if ploidy == 2:
            manifest.update(ploidy=2, heterozygosity_requested=args.heterozygosity,
                            heterozygosity_actual=len(variants) / reference_length,
                            snp_count=len(variants), haplotype_fastas=draft_paths,
                            haplotype_origins={name: {"contig": original, "haplotype": haplotype}
                                               for name, (original, haplotype) in origins.items()},
                            coverage_per_haplotype={str(haplotype): {key: value / reference_length for key, value in hap_bases.items()}
                                                    for haplotype, hap_bases in haplotype_bases.items()},
                            variant_coordinates="0-based; haplotype 1 carries ref, haplotype 2 carries alt")
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        for relative in markers:
            path = stage / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"Self-contained synthetic read fixtures; seed {args.seed}. "
                            "Regenerate with: python workflow/scripts/generate_synthetic_reads.py "
                            f"OUTPUT_DIR --chromosome-lengths {' '.join(map(str, args.chromosome_lengths))} "
                            f"--gc-content {args.gc_content} --seed {args.seed} "
                            f"--heterozygosity {args.heterozygosity}. "
                            "See manifest.json and truth.tsv.gz in the dataset root. "
                            "Hi-C contacts are a toy GATC model, not experimental contact evidence. "
                            "Sequence files are ignored by Git.\n")
        for source in sorted(stage.rglob("*")):
            if source.is_file():
                destination = args.output / source.relative_to(stage)
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.link(source, destination)
    print(json.dumps({"reads": dict(counts), "coverage": manifest["coverage"], "hic_contacts": dict(contacts)}, indent=2))


if __name__ == "__main__":
    main()
