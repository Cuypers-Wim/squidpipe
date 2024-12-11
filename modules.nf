
/*
 * Concatenate FASTQ files per barcode folder
 */

process CONCATENATE_FASTQ {
    label 'process_low'
    conda params.conda_envs.default_env

    input:
    path fqDir
    val meta

    output:
    tuple path("${meta.folder}_allreads.fastq.gz"), val(meta)

    script:
    """
    cat ${fqDir}/${meta.folder}/*.fastq.gz > ${meta.folder}_allreads.fastq.gz
    """
}

/*
 * Run Kraken2 on fastq files within barcode files
 */

process RUNKRAKEN2 {

    label 'process_low'
    conda params.conda_envs.default_env

    publishDir "results/kraken2", pattern: "*_krakenReport.txt", mode: 'copy', overwrite: false
    
    tag "Kraken on $meta.name"

    input:
    tuple path(reads_combined), val(meta)

    output:
    path "*_krakenReport.txt", emit: 'kraken_report'
    tuple path(reads_combined), path("*_ids.txt"), val(meta), optional: true, emit: 'reads_ids_meta'

    script:
    """

    kraken2 \
        --db ${params.databases.kraken_db} \
        --threads ${task.cpus} \
        --gzip-compressed \
        --minimum-hit-groups 1 \
        --output ${meta.name}_krakenReport.txt \
        ${reads_combined}

        awk -F'\\t' -v taxID="${meta.taxid}" '\$1 == "C" && \$3 == taxID {print \$2}' ${meta.name}_krakenReport.txt > ${meta.name}_ids.txt.temp

    # Only move the tmp file to the final location if it's not empty
    if [ -s ${meta.name}_ids.txt.temp ]; then
        mv ${meta.name}_ids.txt.temp ${meta.name}_ids.txt
    else
        echo "${meta.name}_ids.txt.temp is empty. No file will be outputted."
        # We simply don't create the output file, leading to "nothing" being emitted.
    fi
    """
}

/*
 * Extract a subset fastq containing reads of interest
 */

process SUBSET_FASTQ {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    tuple path(reads), path(virus_ids), val(meta)

    output:
    tuple path("${meta.name}.fastq.gz"), val(meta)

    script:
    """
    seqtk subseq ${reads} ${virus_ids} > ${meta.name}.fastq
    pigz -p ${task.cpus} ${meta.name}.fastq
    """

}

/*
 * Retrieve reference genomes based on Taxon IDs
 */

process RETRIEVE_GENOMEREFS {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/references", pattern: "*.fasta", mode: 'copy', overwrite: false
    publishDir "results/references", pattern: "*.log", mode: 'copy', overwrite: false

    input:
    tuple path(subsetted_fastqs), val(meta)

    output:
    tuple path(subsetted_fastqs), path("*.fasta"), val(meta)

    script:
    """
    genome_retriever.py --taxon_id ${meta.taxid} --name ${meta.name} -o .
    
    awk -v add="${meta.taxid}_${meta.name}_" '/^>/ {print ">" add substr(\$0, 2); next} {print}' *.fna > ${meta.name}.fasta
    """
}



/*
 * Deduplicate reference
 */


process DEDUPLICATE_REFERENCE {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    path fasta

    output:
    path 'unique.fasta'

    script:
    """
    seqkit rmdup -n ${fasta} > unique.fasta
    """

}

/*
 * Map reads and perform samtools sort on the resulting bam file
 */

process MAP_READS {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    tuple path(subsetted_fastq), val(meta)
    path fasta

    output:
    tuple path("*.sorted.bam"), val(meta)

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${fasta} ${subsetted_fastq} > mapped.sam
    samtools view -@ ${task.cpus} -b mapped.sam -o mapped.bam
    samtools sort -l 0 -m 1G -O BAM -@ ${task.cpus} mapped.bam -o ${meta.name}_${meta.taxid}.sorted.bam
    """

}

/*
 * Mark duplicates and index the bam files
 */
process BAM_PROCESSING {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    tuple path(bam_sorted), val(meta)

    output:
    tuple path ("*.sorted.dedup.bam"), path ("*.sorted.dedup.bam.bai"), val(meta)

    script:
    """
    samtools markdup -r -@ ${task.cpus} ${bam_sorted} ${meta.name}.sorted.dedup.bam
    samtools index ${meta.name}.sorted.dedup.bam
    """

}

/*
 * Obtain mapping stats on the maping to all genomes
 */

process MAPPING_STATS_ALL {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/mapping", pattern: "*.flagstat.txt", mode: 'copy', overwrite: false
    publishDir "results/mapping", pattern: "*.idxstats.txt", mode: 'copy', overwrite: false
    
    input:
    tuple path (bam_dedup), path (bam_dedup_index), val(meta)

    output:
    path "*.flagstat.txt", emit: 'flags'
    path "*.idxstats.txt", emit: 'stats'
    tuple path(bam_dedup), path(bam_dedup_index), val(meta), emit: 'tuple_bam_meta'

    script:
    """
    samtools flagstat ${bam_dedup} > ${meta.name}.flagstat.txt
    samtools idxstats ${bam_dedup} > ${meta.name}.idxstats.txt
    """
}

/*
 * Extract reads that map to the species of interest
 */
 
process EXTRACT_READIDS {
    
    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/mapping", pattern: "*_filtered.bam", mode: 'copy', overwrite: false
    publishDir "results/mapping", pattern: "*_final_read_ids.txt", mode: 'copy', overwrite: false
    publishDir "results/reads", pattern: "*.fastq.gz", mode: 'copy', overwrite: false

    input:
    tuple path (bam_dedup), path (bam_dedup_index), val(meta)

    output:
    tuple path("*_filtered.bam"), val(meta), emit: 'tuple_bam_meta'
    tuple path("*_final_read_ids.txt"), val(meta), emit: 'tuple_ids_meta'
    tuple path("*.fastq.gz"), val(meta), emit: 'filtered_tp_reads'

    script:
    """
    # extract the correct region name. 

    REGION=\$(samtools view -h ${bam_dedup} | \
    grep "^@SQ" | grep "${meta.taxid}_${meta.name}" | cut -f2 | sed 's/SN://g')

    # use REGION to filter bam wile retaining header
    
    samtools view -b -h -F 256 -F 4 ${bam_dedup} -o ${meta.taxid}_${meta.name}_filtered.bam \$REGION
    samtools view ${meta.taxid}_${meta.name}_filtered.bam | cut -f 1 > ${meta.taxid}_${meta.name}_final_read_ids.txt
    
    # make final fastq with true positives from bam to be pulished 
    samtools fastq -@ ${task.cpus} ${meta.taxid}_${meta.name}_filtered.bam > ${meta.taxid}_${meta.name}.fastq
    pigz -p ${task.cpus} ${meta.taxid}_${meta.name}.fastq
    """
}

/*
 * Run samtools depth on the subsetted bam files
 */

process SAMTOOLS_DEPTH {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/mapping", pattern: "*.depth", mode: 'copy', overwrite: false

    input:
    tuple path(filtered_bam), val(meta)

    output:
    tuple path("*.depth"), val(meta)
    
    script:
    """
    samtools depth -aa ${filtered_bam} | grep "${meta.taxid}_${meta.name}" > ${meta.taxid}_${meta.name}.depth
    """

}

/*
 * Calculate depth of coverage statistics
 */

process DEPTH_OF_COVERAGE {
    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results", pattern: "*.csv", mode: 'copy', overwrite: false

    input:
    path(depthFile)

    output:
    path "coverage_metrics.csv"
    
    script:
    """
    covstats.py ${depthFile} coverage_metrics.csv
    """

}

/*
 * Subset POD5 files based on the target read IDs extracted from the bam file
 */

process SUBSET_POD5 {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/pod5", pattern: "*.filtered.pod5", mode: 'copy', overwrite: false

    input:
    path pod5_reads
    tuple path(virus_ids), val(meta)

    output:
    path("*.filtered.pod5")

    script:
    """
    # Check if there are any subfolders 
    if find ${pod5_reads} -mindepth 1 -type d | read; then
        # Subfolders are present, run this command
        pod5 filter ${pod5_reads}/${meta.folder} --output ${meta.taxid}_${meta.name}.filtered.pod5 --ids ${virus_ids} --missing-ok
    else
        # No subfolders present, just query all files
        pod5 filter ${pod5_reads} --output ${meta.taxid}_${meta.name}.filtered.pod5 --ids ${virus_ids} --missing-ok
    fi

    """

}
