import sys

with open('sample_names.txt') as f:
    SAMPLES = f.read().splitlines()

print('samples are:', SAMPLES, file=sys.stderr)

with open('references.txt') as f2:
    REF = f2.read().splitlines()
    REF_g = [REF[0]]
    REF_r = [REF[1]]

rule all:
    input:
        expand("results/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("results/fastp/{sample}_1P.fq", sample=SAMPLES),
        expand("results/fastp/{sample}_2P.fq", sample=SAMPLES),
	"results/Salmon_index",
        expand("results/Salmon_quant/{sample}.quant", sample=SAMPLES),
        "results/Salmon_quant/all_quant.sf"

rule make_fastqc:
    input:
        R1 = "{sample}_R1_001.fastq.gz",
	R2 = "{sample}_R2_001.fastq.gz",
    params:
        outdir = "results/fastqc",
    priority: 2
    log:
        R1 = "logs/fastqc/{sample}_R1.txt",
        R2 = "logs/fastqc/{sample}_R2.txt",
    output:
        "results/fastqc/{sample}_R1_001_fastqc.html",
        "results/fastqc/{sample}_R1_001_fastqc.zip",
	"results/fastqc/{sample}_R2_001_fastqc.html",
        "results/fastqc/{sample}_R2_001_fastqc.zip",
    run:
        shell("fastqc {input.R1} -o results/fastqc 2> {log.R1}")
	shell("fastqc {input.R2} -o results/fastqc 2> {log.R2}")


rule trim_fastp:
    input:
        R1 = "{sample}_R1_001.fastq.gz",
        R2 = "{sample}_R2_001.fastq.gz",
    params:
        outdir = "results/fastp",
    priority: 2	
    log:
        "logs/fastp/{sample}.txt",
    output:
        outR1 = "results/fastp/{sample}_1P.fq",
        outR2 = "results/fastp/{sample}_2P.fq",
    shell:
        "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.outR1} --out2 {output.outR2} --detect_adapter_for_pe 2> {log}"

rule index_reference:
    input:
        g = REF_g,
        r = REF_r,
    output:
        directory("results/Salmon_index")
    log:
        "logs/salmon/index.txt",
    priority: 3
    run:
        shell("grep '^>' <(gunzip -c {input.g}) | cut -d ' ' -f 1 > decoys.txt")
        shell("sed -i.bak -e 's/>//g' decoys.txt")
        shell("cat {input.r} {input.g} > gentrome.fa.gz")
        shell("salmon index --index {output} --transcripts gentrome.fa.gz -d decoys.txt --gencode 2> {log}")

rule salmon_quant:
    input: 
        R1 = "results/fastp/{sample}_1P.fq",
        R2 = "results/fastp/{sample}_2P.fq",
    log:
        "logs/salmon/{sample}_quant.txt",
    params:
        outdir = "results/Salmon_quant",
    priority: 1
    output: 
        directory("results/Salmon_quant/{sample}.quant"),
    shell:
        "salmon quant -i results/Salmon_index --libType A -1 {input.R1} -2 {input.R2} -o {output} --validateMappings --seqBias --gcBias 2> {log}"

rule salmon_merge:
    input: 
        expand("results/Salmon_quant/{sample}.quant", sample=SAMPLES),
    log:
        "logs/salmon/merge.txt",
    output: 
        "results/Salmon_quant/all_quant.sf"
    priority: 0
    shell:
        "salmon quantmerge --quants {input} -o {output} 2> {log}"
