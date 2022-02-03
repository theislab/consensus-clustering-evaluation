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

    publishDir "$params.outdir/rds_datasets"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("${file.baseName}.Rds"), val(labels)

    script:
    """
    convert_h5ad_Rds.R --out-file "${file.baseName}.Rds" $file
    """
}

process RDS2H5AD {
    conda "envs/zellkonverter.yml"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        tuple val(name), path("${method}.h5ad"), val(labels), val(method)

    script:
    """
    convert_Rds_h5ad.R --out-file "${method}.h5ad" $file
    """
}

process METHOD_SCANPY {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

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

    publishDir "$params.outdir/method_output/${name}"

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
        path("ari_${name}_${method}.tsv")

    script:
    """
    metric_ari.py \\
        --out-file "ari_${name}_${method}.tsv" \\
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
        path("fmi_${name}_${method}.tsv")
    script:
    """
    metric_fmi.py \\
        --out-file "fmi_${name}_${method}.tsv" \\
        --dataset "$name" \\
        --labels $labels \\
        --method $method \\
        $file
    """
}

process COMBINE_METRICS {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir", mode: "copy"

    input:
        path(files)

    output:
        path("metrics.tsv")

    script:
    """
    combine_metrics.py --out-file metrics.tsv $files
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
    metrics_ch = METRIC_ARI.out
        .concat(METRIC_FMI.out)
        .toList()
    COMBINE_METRICS(metrics_ch)
}

/*
========================================================================================
    THE END
========================================================================================
*/
