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

process METHOD_RANDOM {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("random.h5ad"), val(labels), val("random")

    script:
    """
    method_random.py --out-file random.h5ad --labels $labels $file
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

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("seurat.Rds"), val(labels), val("Seurat")

    script:
    """
    method_seurat.R --out-file seurat.Rds $file
    """
}

process METHOD_SIMLR {
    conda "envs/simlr.yml"

    label "process_low"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("simlr.Rds"), val(labels), val("SIMLR")

    script:
    """
    method_simlr.R --out-file simlr.Rds --ncpus ${task.cpus} $file
    """
}

process METHOD_SC3 {
    conda "envs/sc3.yml"

    label "process_low"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("sc3.Rds"), val(labels), val("SC3")

    script:
    """
    method_sc3.R --out-file sc3.Rds --ncpus ${task.cpus} $file
    """
}

process RUN_CONSTCLUST {
    conda "envs/constclust.yml"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("constclust.h5ad"), val(labels)

    script:
    """
    run_constclust.py --out-file constclust.h5ad $file
    """
}

process METHOD_MRCC {
    conda "envs/mrcc.yml"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("mrcc.h5ad"), val(labels), val("MRCC")

    script:
    """
    method_mrcc.py --out-file mrcc.h5ad $file
    """
}

process METHOD_COLA {
    conda "envs/cola.yml"

    publishDir "$params.outdir/method_output/${name}", mode: "copy"

    input:
        tuple val(name), path(file), val(labels)

    output:
        tuple val(name), path("cola.Rds"), val(labels), val("cola")

    script:
    """
    method_cola.R --out-file cola.Rds --labels $labels $file
    """
}

process MATCH_CLUSTERS {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/method_output/${name}"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        tuple val(name), path("${method}_matched.h5ad"), val(labels), val(method)

    script:
    """
    match_clusters.py \\
        --clusters "Cluster" \\
        --labels "$labels" \\
        --out-file "${method}_matched.h5ad" \\
        $file
    """
}

process H5AD2RDS_MATCHED {
    conda "envs/zellkonverter.yml"

    publishDir "$params.outdir/method_output/${name}"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        tuple val(name), path("${file.baseName}.Rds"), val(labels), val(method)

    script:
    """
    convert_h5ad_Rds.R --out-file "${file.baseName}.Rds" $file
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

process METRIC_AMI {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        path("ami_${name}_${method}.tsv")

    script:
    """
    metric_ami.py \\
        --out-file "ami_${name}_${method}.tsv" \\
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

process METRIC_COMPLETENESS {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        path("completeness_${name}_${method}.tsv")
    script:
    """
    metric_completeness.py \\
        --out-file "completeness_${name}_${method}.tsv" \\
        --dataset "$name" \\
        --labels $labels \\
        --method $method \\
        $file
    """
}

process METRIC_HOMOGENEITY {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        path("homogeneity_${name}_${method}.tsv")
    script:
    """
    metric_homogeneity.py \\
        --out-file "homogeneity_${name}_${method}.tsv" \\
        --dataset "$name" \\
        --labels $labels \\
        --method $method \\
        $file
    """
}

process METRIC_F1 {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        path("F1_${name}_${method}.tsv")
    script:
    """
    metric_F1.py \\
        --out-file "F1_${name}_${method}.tsv" \\
        --dataset "$name" \\
        --labels $labels \\
        --method $method \\
        $file
    """
}

process METRIC_MCC {
    conda "envs/sklearn.yml"

    input:
        tuple val(name), path(file), val(labels), val(method)

    output:
        path("mcc_${name}_${method}.tsv")
    script:
    """
    metric_mcc.py \\
        --out-file "mcc_${name}_${method}.tsv" \\
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

process PLOT_METRICS {
    conda "envs/tidyverse.yml"

    publishDir "$params.outdir", mode: "copy"

    input:
        path(file)

    output:
        path("metrics-plots.pdf")

    script:
    """
    plot_metrics.R --out-file metrics-plots.pdf $file
    """
}

workflow CCEVAL {
    PROFILE(datasets_ch)
    H5AD2RDS(datasets_ch)
    METHOD_RANDOM(datasets_ch)
    METHOD_SCANPY(datasets_ch)
    METHOD_SEURAT(H5AD2RDS.out)
    METHOD_SIMLR(H5AD2RDS.out)
    METHOD_COLA(H5AD2RDS.out)
    METHOD_SC3(H5AD2RDS.out)
    RUN_CONSTCLUST(datasets_ch)
    METHOD_MRCC(RUN_CONSTCLUST.out)
    rds_ch = METHOD_SEURAT.out
        .concat(METHOD_SIMLR.out)
        .concat(METHOD_SC3.out)
        .concat(METHOD_COLA.out)
    RDS2H5AD(rds_ch)
    output_ch = METHOD_RANDOM.out
        .concat(METHOD_SCANPY.out)
        .concat(METHOD_MRCC.out)
        .concat(RDS2H5AD.out)
    MATCH_CLUSTERS(output_ch)
    H5AD2RDS_MATCHED(MATCH_CLUSTERS.out)
    METRIC_ARI(MATCH_CLUSTERS.out)
    METRIC_AMI(MATCH_CLUSTERS.out)
    METRIC_FMI(MATCH_CLUSTERS.out)
    METRIC_COMPLETENESS(MATCH_CLUSTERS.out)
    METRIC_HOMOGENEITY(MATCH_CLUSTERS.out)
    METRIC_F1(MATCH_CLUSTERS.out)
    METRIC_MCC(MATCH_CLUSTERS.out)
    metrics_ch = METRIC_ARI.out
        .concat(METRIC_AMI.out)
        .concat(METRIC_FMI.out)
        .concat(METRIC_COMPLETENESS.out)
        .concat(METRIC_HOMOGENEITY.out)
        .concat(METRIC_F1.out)
        .concat(METRIC_MCC.out)
        .toList()
    COMBINE_METRICS(metrics_ch)
    PLOT_METRICS(COMBINE_METRICS.out)
}

/*
========================================================================================
    THE END
========================================================================================
*/
