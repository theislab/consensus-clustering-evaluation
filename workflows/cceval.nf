/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

datasets_ch = Channel
    .fromList(params.input).
    map { input ->
        tuple(
            input.name,
            file(params.datasetsdir + "/" + input.file),
            input.labels
        )
    }

process PROFILE {
    conda "envs/scanpy.yml"

    input:
        tuple val(name), path(file), val(labels)

    script:
    """
    profile_dataset.py --name "$name" --labels "$labels" $file
    """
}

process METHOD_SCANPY {
    conda "envs/scanpy.yml"

    input:
        tuple val(name), path(file), val(labels)

    output:
        path("scanpy.h5ad")

    script:
    """
    method_scanpy.py --out-file scanpy.h5ad $file
    """
}

workflow CCEVAL {
    PROFILE(datasets_ch)
    METHOD_SCANPY(datasets_ch)
}

/*
========================================================================================
    THE END
========================================================================================
*/
