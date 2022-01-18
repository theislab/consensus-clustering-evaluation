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

process H5AD2RDS {
    conda "envs/zellkonverter.yml"

    input:
        tuple val(name), path(file), val(labels)

    output:
        path("dataset.Rds")

    script:
    """
    convert_h5ad_Rds.R --out-file dataset.Rds $file
    """
}

process RDS2H5AD {
    conda "envs/zellkonverter.yml"

    input:
        path(file)

    output:
        path("dataset.h5ad")

    script:
    """
    convert_Rds_h5ad.R --out-file dataset.h5ad $file
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

process METHOD_SEURAT {
    conda "envs/seurat.yml"

    input:
        path(file)

    output:
        path("seurat.Rds")

    script:
    """
    method_seurat.R --out-file seurat.Rds $file
    """
}

workflow CCEVAL {
    PROFILE(datasets_ch)
    H5AD2RDS(datasets_ch)
    METHOD_SCANPY(datasets_ch)
    METHOD_SEURAT(H5AD2RDS.out)
    rds_ch = METHOD_SEURAT.out
    RDS2H5AD(rds_ch)
}

/*
========================================================================================
    THE END
========================================================================================
*/
