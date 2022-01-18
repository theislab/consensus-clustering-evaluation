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
        tuple val(name), path("dataset.Rds"), val(labels)

    script:
    """
    convert_h5ad_Rds.R --out-file dataset.Rds $file
    """
}

process RDS2H5AD {
    conda "envs/zellkonverter.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        tuple val(name), path("dataset.h5ad"), val(labels), val(method)

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
        tuple val(name), path("scanpy.h5ad"), val(labels), val("scanpy")

    script:
    """
    method_scanpy.py --out-file scanpy.h5ad $file
    """
}

process METHOD_SEURAT {
    conda "envs/seurat.yml"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("seurat.Rds"), val(labels), val("Seurat")

    script:
    """
    method_seurat.R --out-file seurat.Rds $file
    """
}

process METRIC_ARI {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        tuple val(name), path("ari.tsv"), val(labels), val(method)

    script:
    """
    metric_ari.py \\
        --out-file ari.tsv \\
        --dataset "$name" \\
        --labels $labels \\
        --method $method \\
        $file
    """
}

process METRIC_FMI {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        tuple val(name), path("fmi.tsv"), val(labels), val(method)

    script:
    """
    metric_fmi.py \\
        --out-file fmi.tsv \\
        --dataset "$name" \\
        --labels $labels \\
        --method $method \\
        $file
    """
}

workflow CCEVAL {
    PROFILE(datasets_ch)
    H5AD2RDS(datasets_ch)
    METHOD_SCANPY(datasets_ch)
    METHOD_SEURAT(H5AD2RDS.out)
    rds_ch = METHOD_SEURAT.out
    RDS2H5AD(rds_ch)
    output_ch = METHOD_SCANPY.out.concat(RDS2H5AD.out)
    METRIC_ARI(output_ch)
    METRIC_FMI(output_ch)
}

/*
========================================================================================
    THE END
========================================================================================
*/
