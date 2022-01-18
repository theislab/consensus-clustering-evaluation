#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/cceval
========================================================================================
    Github : https://github.com/nf-core/cceval
    Website: https://nf-co.re/cceval
    Slack  : https://nfcore.slack.com/channels/cceval
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { CCEVAL } from './workflows/cceval'

//
// WORKFLOW: Run main nf-core/cceval analysis pipeline
//
workflow NFCORE_CCEVAL {
    CCEVAL ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_CCEVAL ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
