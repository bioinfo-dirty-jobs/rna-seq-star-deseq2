def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


rule align:
    input:
        sample=get_trimmed,
        dt=config['ref']['index']+"chrLength.txt",
    output:
        # see STAR manual for additional output files
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 24
    wrapper:
        "0.19.4/bio/star/align"



rule star_index:
    input:
        config['ref']['transcriptome_fasta']
    output:
        dt=config['ref']['index']+"chrLength.txt"
    threads: 1
    params:
        gtf = config['ref']['annotation'],
        genome=config['ref']['index'],
    conda:
    	"envs/star.yaml"
    shell:
          'STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.genome} --genomeFastaFiles {input}  --sjdbGTFfile {params.gtf} --sjdbOverhang 49'
