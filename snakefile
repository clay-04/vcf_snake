import os
from glob import glob
#sample configuration
SAMPLES = sorted(set(
    os.path.basename(f).replace('_L007_R1_001.fastq.gz', '')
    for f in glob("/rhome/clay007/bigdata/BWA/raw_data/WAA/*_L007_R1_001.fastq.gz")
))

#reference genome
ref = "/rhome/clay007/bigdata/BWA/reference_data/Eriosoma_lanigerum_v1.0.fa"

#define readgroups
def get_read_group(wildcards):
    return f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:{wildcards.sample}_lib1\\tPU:unit1"


#CPU and parallel config
TOTAL_CPUS = 100
TOTAL_MEM = 500000 #500gb in mb

#parallel
SAMPLES_IN_PARALLEL = 4 #align/sort/markdups
HC_SAMPLES_IN_PARALLEL = 2 #haplotypecaller

#memory
ALIGN_THREADS = 20
MARKDUPS_THREADS = 10
HC_THREADS = 20


rule all:
    input:
        expand("/rhome/clay007/bigdata/BWA/waa/gatk/{sample}_var.g.vcf.gz", sample=SAMPLES),
        expand("/rhome/clay007/bigdata/BWA/waa/gatk/{sample}_var.g.vcf.gz.tbi", sample=SAMPLES),
        "/rhome/clay007/bigdata/BWA/waa/qc/master_summary.csv"
rule fastp_trim:
    input:
        r1 = "/rhome/clay007/bigdata/BWA/raw_data/WAA/{sample}_L007_R1_001.fastq.gz",
        r2 = "/rhome/clay007/bigdata/BWA/raw_data/WAA/{sample}_L007_R2_001.fastq.gz"
    output:
        r1 = "/rhome/clay007/bigdata/BWA/waa/trim/{sample}_R1_trimmed.fastq.gz",
        r2 = "/rhome/clay007/bigdata/BWA/waa/trim/{sample}_R2_trimmed.fastq.gz",
        html = "/rhome/clay007/bigdata/BWA/waa/trim/{sample}_fastp.html",
        json = "/rhome/clay007/bigdata/BWA/waa/trim/{sample}_fastp.json"
    threads: 8
    resources: mem_mb=32000
    shell:
        """
        module load fastp/0.23.2
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.html} -j {output.json} \
              --detect_adapter_for_pe \
              --trim_poly_g \
              -w {threads}
        """

rule bwa_map:
    input:
        r1 = "/rhome/clay007/bigdata/BWA/waa/trim/{sample}_R1_trimmed.fastq.gz",
        r2 = "/rhome/clay007/bigdata/BWA/waa/trim/{sample}_R2_trimmed.fastq.gz",
        ref = ref
    output:
        "/rhome/clay007/bigdata/BWA/waa/sort/{sample}_sort.bam"
    params:
        rg = get_read_group
    threads:
        ALIGN_THREADS
    resources: mem_mb=160000
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.19.2
        bwa mem -M -t {threads} -R "{params.rg}" {input.ref} {input.r1} {input.r2} | samtools sort -o {output}
        """

rule mark_dups:
    input:
        "/rhome/clay007/bigdata/BWA/waa/sort/{sample}_sort.bam"
    output:
        bam = "/rhome/clay007/bigdata/BWA/waa/dups/{sample}_dup.bam",
        bai = "/rhome/clay007/bigdata/BWA/waa/dups/{sample}_dup.bam.bai"
    threads:
        MARKDUPS_THREADS
    resources: mem_mb=120000
    shell:
        """
        module unload java
        module load java/17.0.2
        module load gatk/4.4.0.0
        module load workspace/scratch
        gatk MarkDuplicatesSpark -I {input} -O {output.bam} --conf "spark.local.dir=$SCRATCH" \
        """

rule flagstat_qc:
    input:
        "/rhome/clay007/bigdata/BWA/waa/dups/{sample}_dup.bam"
    output:
        "/rhome/clay007/bigdata/BWA/waa/qc/{sample}_flagstat.txt"
    shell:
        """
        module load samtools/1.19.2
        samtools flagstat {input} > {output}
        """

rule coverage_qc:
    input:
        "/rhome/clay007/bigdata/BWA/waa/dups/{sample}_dup.bam"
    output:
        "/rhome/clay007/bigdata/BWA/waa/qc/{sample}_coverage.txt"
    shell:
        """
        module load samtools/1.19.2
        samtools coverage {input} > {output}
        """

rule haplotype_caller:
    input:
        dup="/rhome/clay007/bigdata/BWA/waa/dups/{sample}_dup.bam",
        bai="/rhome/clay007/bigdata/BWA/waa/dups/{sample}_dup.bam.bai",
        ref=ref
    output:
        gvcf="/rhome/clay007/bigdata/BWA/waa/gatk/{sample}_var.g.vcf.gz",
        tbi="/rhome/clay007/bigdata/BWA/waa/gatk/{sample}_var.g.vcf.gz.tbi"
    threads:
        HC_THREADS
    resources: mem_mb=220000
    shell:
        """
        module load gatk/4.4.0.0
        gatk --java-options "-Xmx120G -XX:ParallelGCThreads={threads}" HaplotypeCaller \
        -R {input.ref} -I {input.dup} \
        --sample-name {wildcards.sample} \
        -O {output.gvcf} -ERC GVCF
        """

rule bcftools_stats:
    input:
        "/rhome/clay007/bigdata/BWA/waa/gatk/{sample}_var.g.vcf.gz"
    output:
        "/rhome/clay007/bigdata/BWA/waa/qc/{sample}.vcf.stats"
    shell:
        """
        module load bcftools
        bcftools stats {input} > {output}
        """

rule collect_summary_stats:
    input:
        fastp_json = lambda wildcards: glob("/rhome/clay007/bigdata/BWA/waa/trim/*_fastp.json"),
        flagstat = lambda wildcards: glob("/rhome/clay007/bigdata/BWA/waa/qc/*_flagstat.txt"),
        coverage = lambda wildcards: glob("/rhome/clay007/bigdata/BWA/waa/qc/*_coverage.txt"),
        vcfstats = lambda wildcards: glob("/rhome/clay007/bigdata/BWA/waa/qc/*.vcf.stats")
    output:
        "/rhome/clay007/bigdata/BWA/waa/qc/master_summary.csv"
    run:
        import json
        import pandas as pd
        import re

        # Initialize empty list for all sample summaries
        all_stats = []

        # Helper to extract sample name from filepath
        def sample_name(f):
            return re.search(r"/([^/]+?)(?:_fastp.json|_flagstat.txt|_coverage.txt|.vcf.stats)$", f).group(1)

        # Process fastp JSON
        for f in input.fastp_json:
            with open(f) as jf:
                data = json.load(jf)
            sample = sample_name(f)
            stats = {
                "sample": sample,
                "fastp_total_reads": int(data['summary']['before_filtering']['total_reads']),
                "fastp_q30_rate": float(data['summary']['before_filtering']['q30_rate'])
            }
            all_stats.append(stats)

        # Process flagstat and coverage
        for f in input.flagstat:
            sample = sample_name(f)
            with open(f) as ff:
                lines = ff.readlines()
            total_reads = int(lines[0].split()[0])
            # Find existing sample dict and update
            for s in all_stats:
                if s["sample"] == sample:
                    s["flagstat_total_reads"] = total_reads

        for f in input.coverage:
            sample = sample_name(f)
            with open(f) as cf:
                lines = cf.readlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                avg_coverage = float(line.split()[6])
                break  # samtools coverage first line, 7th column = mean coverage
            for s in all_stats:
                if s["sample"] == sample:
                    s["avg_coverage"] = avg_coverage

        # Process VCF stats
        for f in input.vcfstats:
            sample = sample_name(f)
            qual_bins = []
            ts_tv_bins = []
            snp_count = 0

            with open(f) as vf:
                for l in vf:
                    if l.startswith("SN"):
                        if "number of SNPs:" in l:
                            snp_count = int(l.strip().split(":")[1].strip())
                    elif l.startswith("QUAL"):
                        parts = l.strip().split("\t")
                        if len(parts) < 7:
                            continue
                        qual_bin = float(parts[2])
                        num_snps = int(parts[3])
                        ts = int(parts[5])
                        tv = int(parts[6])

                        # Collect for mean qual calculation
                        qual_bins.append((qual_bin,num_snps))

                        # Ts/Tv per bin
                        ts_tv_ratio = ts / tv if tv != 0 else None
                        ts_tv_bins.append((qual_bin, ts_tv_ratio, num_snps))

        # Compute mean quality weighted by SNP counts
            total_snps = sum(x[1] for x in qual_bins)
            mean_qual = sum(q*num for q, num in qual_bins) / total_snps if total_snps != 0 else 0

        #update the all_stats dict for this sample
            for s in all_stats:
                if s["sample"] == sample:
                    s["snp_count"] = snp_count
                    s["mean_qual"] = mean_qual
                    s["qual_bins"] = qual_bins
                    s["ts_tv_bins"] = ts_tv_bins


        # Convert to DataFrame and write CSV
        df = pd.DataFrame(all_stats)
        df.to_csv(output[0], index=False)

