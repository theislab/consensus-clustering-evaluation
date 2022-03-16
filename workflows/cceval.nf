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

mrcc_ch = Channel
    .fromList(params.mrcc).
    map { mrcc ->
        tuple(
            mrcc.name,
            mrcc.graph_type,
            mrcc.community_type,
            mrcc.outlier_type,
            mrcc.outlier_thresh,
            mrcc.merge_thresh
        )
    }

process PREP {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/prepped_datasets"

    input:
        tuple val(dataset), path(file), val(labels)

    output:
        tuple val(dataset), path("${dataset}.h5ad")

    script:
    """
    prep_dataset.py --name "$dataset" --labels "$labels" --out-file "${dataset}.h5ad" $file
    """
}

process PREP_RDS {
    conda "envs/zellkonverter.yml"

    publishDir "$params.outdir/prepped_datasets"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("${file.baseName}.Rds")

    script:
    """
    convert_h5ad_Rds.R --out-file "${file.baseName}.Rds" $file
    """
}

process METHOD_RANDOM {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("random.tsv"), val("random")

    script:
    """
    method_random.py --out-file random.tsv $file
    """
}

process METHOD_SCANPY {
    conda "envs/scanpy.yml"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("scanpy.tsv"), val("scanpy")

    script:
    """
    method_scanpy.py --out-file scanpy.tsv $file
    """
}

process METHOD_SEURAT {
    conda "envs/seurat.yml"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("seurat.tsv"), val("Seurat")

    script:
    """
    method_seurat.R --out-file seurat.tsv $file
    """
}

process METHOD_SIMLR {
    conda "envs/simlr.yml"

    label "process_low"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("simlr.tsv"), val("SIMLR")

    script:
    """
    method_simlr.R --out-file simlr.tsv --ncpus ${task.cpus} $file
    """
}

process METHOD_SC3 {
    conda "envs/sc3.yml"

    label "process_medium"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("sc3.tsv"), val("SC3")

    script:
    """
    method_sc3.R --out-file sc3.tsv --ncpus ${task.cpus} $file
    """
}

process METHOD_COLA {
    conda "envs/cola.yml"

    label "process_high"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("cola.tsv"), val("cola")

    script:
    """
    method_cola.R --out-file cola.tsv --ncpus ${task.cpus} $file
    """
}

process RUN_CONSTCLUST {
    conda "envs/constclust.yml"

    publishDir "$params.outdir/method_output/${dataset}"

    input:
        tuple val(dataset), path(file)

    output:
        tuple val(dataset), path("constclust.h5ad")

    script:
    """
    run_constclust.py --out-file constclust.h5ad $file
    """
}

process METHOD_MRCC {
    conda "envs/mrcc.yml"

    input:
        tuple val(dataset), path(file), val(name), val(graph_type), val(community_type),
        val(outlier_type), val(outlier_thresh), val(merge_thresh)

    output:
        tuple val(dataset), path("${name}.tsv"), val(name)

    script:
    """
    method_mrcc.py \\
        --out-file ${name}.tsv \\
        --graph-type ${graph_type} \\
        --community-type ${community_type} \\
        --outlier-type ${outlier_type} \\
        --outlier-thresh ${outlier_thresh} \\
        --merge-thresh ${merge_thresh} \\
        $file
    """
}

process MATCH_CLUSTERS {
    conda "envs/sklearn.yml"

    publishDir "$params.outdir/method_output/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        tuple val(dataset), path("${method}_clusters.tsv"), val(method)

    script:
    """
    match_clusters.py --out-file "${method}_clusters.tsv" $file
    """
}

process METRIC_ARI {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("ari_${dataset}_${method}.tsv")

    script:
    """
    metric_ari.py \\
        --out-file "ari_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_AMI {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("ami_${dataset}_${method}.tsv")

    script:
    """
    metric_ami.py \\
        --out-file "ami_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_FMI {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("fmi_${dataset}_${method}.tsv")

    script:
    """
    metric_fmi.py \\
        --out-file "fmi_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_COMPLETENESS {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("completeness_${dataset}_${method}.tsv")

    script:
    """
    metric_completeness.py \\
        --out-file "completeness_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_HOMOGENEITY {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("homogeneity_${dataset}_${method}.tsv")

    script:
    """
    metric_homogeneity.py \\
        --out-file "homogeneity_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_F1 {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("F1_${dataset}_${method}.tsv")
    script:
    """
    metric_F1.py \\
        --out-file "F1_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_MCC {
    conda "envs/sklearn.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("mcc_${dataset}_${method}.tsv")
    script:
    """
    metric_mcc.py \\
        --out-file "mcc_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
        $file
    """
}

process METRIC_ECS {
    conda "envs/clustassess.yml"

    input:
        tuple val(dataset), path(file), val(method)

    output:
        path("ecs_${dataset}_${method}.tsv")
    script:
    """
    metric_ecs.R \\
        --out-file "ecs_${dataset}_${method}.tsv" \\
        --dataset "$dataset" \\
        --method "$method" \\
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

process PLOT_UMAP {
    conda "envs/scanpy.yml"

    publishDir "$params.outdir/method_output/${dataset[0]}/umaps", mode: "copy"

    input:
        tuple val(dataset), val(clusters)

    output:
        path("${clusters[2]}.pdf")

    script:
    """
    plot_umap.py \\
        --out-file "${clusters[2]}.pdf" \\
        --clusters "${clusters[1]}" \\
        "${dataset[1]}"
    """
}

workflow CCEVAL {
    PREP(datasets_ch)
    PREP_RDS(PREP.out)
    METHOD_RANDOM(PREP.out)
    METHOD_SCANPY(PREP.out)
    METHOD_SEURAT(PREP_RDS.out)
    METHOD_SIMLR(PREP_RDS.out)
    METHOD_COLA(PREP_RDS.out)
    METHOD_SC3(PREP_RDS.out)
    RUN_CONSTCLUST(PREP.out)
    METHOD_MRCC(RUN_CONSTCLUST.out.combine(mrcc_ch))
    output_ch = METHOD_RANDOM.out
        .concat(METHOD_SCANPY.out)
        .concat(METHOD_MRCC.out)
        .concat(METHOD_SEURAT.out)
        .concat(METHOD_SIMLR.out)
        .concat(METHOD_SC3.out)
        .concat(METHOD_COLA.out)
    MATCH_CLUSTERS(output_ch)
    METRIC_ARI(MATCH_CLUSTERS.out)
    METRIC_AMI(MATCH_CLUSTERS.out)
    METRIC_FMI(MATCH_CLUSTERS.out)
    METRIC_COMPLETENESS(MATCH_CLUSTERS.out)
    METRIC_HOMOGENEITY(MATCH_CLUSTERS.out)
    METRIC_F1(MATCH_CLUSTERS.out)
    METRIC_MCC(MATCH_CLUSTERS.out)
    METRIC_ECS(MATCH_CLUSTERS.out)
    metrics_ch = METRIC_ARI.out
        .concat(METRIC_AMI.out)
        .concat(METRIC_FMI.out)
        .concat(METRIC_COMPLETENESS.out)
        .concat(METRIC_HOMOGENEITY.out)
        .concat(METRIC_F1.out)
        .concat(METRIC_MCC.out)
        .concat(METRIC_ECS.out)
        .toList()
    COMBINE_METRICS(metrics_ch)
    PLOT_METRICS(COMBINE_METRICS.out)
    PLOT_UMAP(PREP.out.cross(MATCH_CLUSTERS.out))
}

/*
========================================================================================
    THE END
========================================================================================
*/
