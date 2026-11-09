nextflow.enable.dsl=2

// Execute DRAM2's original process, without gene calling or unrelated searches.
include { RUNDBCAN_CAZYMEANNOTATION } from './upstream/module/main'

workflow {
    RUNDBCAN_CAZYMEANNOTATION(
        Channel.of(tuple([id: 'genome'], file(params.fasta, checkIfExists: true))),
        file(params.db, checkIfExists: true)
    )
}
