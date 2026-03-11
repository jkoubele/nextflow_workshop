
process FastQC {

    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(sample), path(read1), path(read2)

    tag "$sample"

    output:
    path("*.zip"), emit: fastqc_reports
    path("*.html")

    publishDir "./results/FastQC/${sample}", mode: 'copy'

    script:
    """
    fastqc $read1 $read2 --outdir .
    """
}

process MultiQC {

    container "multiqc/multiqc:v1.30"

    input:
    path zip_files

    output:
    path("*")


    publishDir "./results/MultiQC", mode: 'copy'

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

workflow{

    samplesheet = './data/samplesheet.csv'
    fastq_dir = './data/FASTQ'
    samples = Channel
            .fromPath(samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def sample      = row.sample
                def fq1         = file("${fastq_dir}/${row.fastq_1}")
                def fq2         = file("${fastq_dir}/${row.fastq_2}")

                if (!fq1.exists()) error "FASTQ file not found: ${fq1}"
                if (!fq2.exists()) error "FASTQ file not found: ${fq2}"

                tuple(sample, fq1, fq2)
            }

        samples | view
        def fastqc_out = FastQC(samples)
        def fastqc_out_aggregated = fastqc_out.fastqc_reports.collect()
        MultiQC(fastqc_out_aggregated)
}
