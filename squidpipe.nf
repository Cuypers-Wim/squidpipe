#!/usr/bin/env nextflow


/*

    S Q U I D P I P E - N F 
    ===================================

    ABOUT
    -----
    SQUIDPIPE-NF is a Nextflow pipeline designed for processing nanopore sequencing data, 
    focusing on extracting raw signals and reads from specific taxa in barcoded runs. 
    The pipeline utilizes Kraken2 for taxonomic classification, extracts reads corresponding 
    to specified taxa, and performs subsequent mapping and depth-of-coverage analysis. 
    The pipeline is modular and scalable, supporting the integration of custom reference genomes 
    and POD5 file processing.


*/


nextflow.enable.dsl=2

// include 

include { CONCATENATE_FASTQ     } from './modules.nf'
include { RUNKRAKEN2            } from './modules.nf'
include { SUBSET_FASTQ          } from './modules.nf'
include { RETRIEVE_GENOMEREFS   } from './modules.nf'
include { DEDUPLICATE_REFERENCE } from './modules.nf'
include { MAP_READS             } from './modules.nf'
include { BAM_PROCESSING        } from './modules.nf'
include { MAPPING_STATS_ALL     } from './modules.nf'
include { EXTRACT_READIDS       } from './modules.nf'
include { SAMTOOLS_DEPTH        } from './modules.nf'
include { DEPTH_OF_COVERAGE     } from './modules.nf'
include { SUBSET_POD5           } from './modules.nf'

// workflow

workflow {
    if( params.csvMeta )
        ch_csv_lines = Channel.fromPath(params.csv_file)
                        .splitCsv(header: true, sep: ";")
                        .map { row ->
                            meta = [
                                folder: "${row.filename}",
                                name: "${row.species_name}",
                                taxid: "${row.pathogen_taxid}",
                                year: "${row.year_of_isolation}",
                                country: "${row.country_of_isolation}",
                                ori: "${row.geographic_origin}",
                                strain: "${row.strain_lineage}",
                                source: "${row.source_id}",
                                host: "${row.host_ncbi_taxid}",
                                labid: "${row.internal_lab_id}",
                                diagMethod: "${row.diagnostic_method_id}",
                                remarks: "${row.remarks}"
                            ]
                        }
    else
        ch_csv_lines = Channel.fromPath(params.csv_file)
                        .splitCsv(header: false)
                        .map { row ->
                            meta = [
                                folder: row[0],
                                name: row[1],
                                taxid: row[2]
                            ]
                        }
                        // produces:
                            // val(meta)
                            // meta.folder
                            // meta.name
                            // meta.taxid

    concatenate_fastq_results = CONCATENATE_FASTQ(params.fastq_dir, ch_csv_lines)
    kraken_results = RUNKRAKEN2(concatenate_fastq_results)

    fastq_subset_results = SUBSET_FASTQ(kraken_results.reads_ids_meta)
    // subset fastq emits a tuple of subsetted reads, and virus ids

    // retrieve genome references
    
    ch_retrieved_refs = RETRIEVE_GENOMEREFS(fastq_subset_results)

    // Split the results into two channels: fna_files and fastq_meta

    ch_mapping_input = ch_retrieved_refs.multiMap{ it ->
        fna_files: it[1]
        reads_meta: it[0, 2]
    }

    // combine all downloaded fna files with the user-provided fna file (typically the "decoy" fasta(s) to avoid false-positives)
    ch_extraRef = channel.fromPath(params.references)

    ch_all_refs = ch_mapping_input.fna_files.concat(ch_extraRef)
    
    ch_refGenome = ch_all_refs.collectFile(name: 'combined_reference.fna')

    // deduplicate and  to a value channel convert ch_refGenome to a value channel

    ch_refGenome_value = DEDUPLICATE_REFERENCE(ch_refGenome).first()

    // map reads per sample to the conbined reference

    mapping_output = MAP_READS(ch_mapping_input.reads_meta, ch_refGenome_value)

    bamProcessing_output = BAM_PROCESSING(mapping_output)

    // todo, how to get correct read ids? 
    // todo, how to get mapsats per correct reference?
    // todo - compare to kraken2 outputs
    
    // extracts the readids of those reads that map to the correct reference genome
    
    mapping_stats_output = MAPPING_STATS_ALL(bamProcessing_output)

    extract_readids_output = EXTRACT_READIDS(mapping_stats_output.tuple_bam_meta)

    // we need only the filtered bam and meta

    samtools_depth_output = SAMTOOLS_DEPTH(extract_readids_output.tuple_bam_meta)

    ch_depthFiles = samtools_depth_output.map{ it -> it[0] }
    ch_concat_deptFile = ch_depthFiles.collectFile(name: 'combined.depth')
    
    depth_of_coverage_output = DEPTH_OF_COVERAGE(ch_concat_deptFile)

    if (params.pod5_split) {
    split_pod5_results = SUBSET_POD5(params.pod5_dir, extract_readids_output.tuple_ids_meta)

    //split_pod5_results.collectFile(name: 'final_metadata.csv', newLine: true, storeDir: 'results/')
    
    //split_pod5_results
    //.map { file, meta -> "${file.getName()},${meta}" }  // Format as CSV
    //.collectFile(name: 'final_metadata.csv', newLine: true, storeDir: 'results/')

    split_pod5_results
    .map { file, meta -> 
        def headers = "filename," + meta.keySet().join(",")  // Add "filename" as the first header
        def values = "${file.getName()}," + meta.values().join(",")
        return [headers, values]
    }
    .flatten()
    .unique()
    .collectFile(name: 'final_metadata.csv', newLine: true, storeDir: 'results/')


    } 

}       

// ToDo
//
// final metadatacsv file should consist of pod5 filename and all metadata
//
// Switch to nf-core modules - better for portability
// docker - apptainer
// Make a nice schema for visual depiction of the pipeline, or print the DAG?
