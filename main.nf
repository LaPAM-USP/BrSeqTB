nextflow.enable.dsl = 2

/*
 * ============================================================
 * PARAMETERS
 * ============================================================
 */

params.run                = null
params.add_kaiju_manually = false
params.input_table        = "input/input_table.xlsx"
params.reads_dir          = "reads"
params.demo = false
params.module = null  // Example: bwa, trimmomatic, cohort, etc.
params.exclude = null

/*
 * ============================================================
 * BLOCK - INIT CHAIN (SEQUENTIAL + CACHED)
 * ============================================================
 */

process KAIJU_DB {

    tag "kaiju-db"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    bash bin/kaijudb.sh ${params.add_kaiju_manually}
    """
}


process OMS_CATALOG {

    tag "oms-catalog"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    python bin/omsCatalog.py
    """
}


process BWA_REF {

    tag "bwa-ref"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    bash bin/bwaref.sh
    """
}


process GATK_DICT {

    tag "gatk-dict"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    bash bin/gatkdict.sh
    """
}


process SNPEFF_DB {

    tag "snpeff-db"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    bash bin/snpeffdb.sh
    """
}


process MAKE_MANIFEST_VALIDATE {

    tag "manifest"

    input:
        val token
        path input_table
        path reads_dir

    output:
        path "manifest.tsv"

    publishDir "${projectDir}", mode: 'copy'

    script:
    """
    python ${projectDir}/bin/make_manifest_validate.py \
        --xlsx ${input_table} \
        --reads ${reads_dir} \
        --out manifest.tsv
    """
}

/*
 * ============================================================
 * BLOCK 1 — PER BIOSAMPLE
 * ============================================================
 */

process FASTQC {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/fastqc.sh ${biosample} ${params.reads_dir}
    """
}


process TRIMMOMATIC {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/trimmomatic.sh ${biosample}
    """
}


process KAIJU {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/kaiju.sh ${biosample}
    """
}


process BWA {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/bwa.sh ${biosample}
    """
}


process DELLY {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/delly.sh ${biosample}
    """
}


process LOFREQ {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/lofreq.sh ${biosample}
    """
}


process GATK_GVCF {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/gatkGvcf.sh ${biosample}
    """
}


process GATK_VCF {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/gatkVcf.sh ${biosample}
    """
}


process NORM {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/norm.sh ${biosample}
    """
}

process TBDR_RCOV {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    python bin/tbdrRCov.py ${biosample}
    """
}

process LINEAGE {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    python bin/lineage.py ${biosample}
    """
}

process NTM_FILTER {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/ntmFilter.sh ${biosample}
    """
}

process SNPEFF {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    bash bin/snpeff.sh ${biosample}
    """
}

/*
 * ============================================================
 * BLOCO 2 — COHORT (GLOBAL)
 * Inicia apenas após sincronização global do BLOCO 1
 * ============================================================
 */

process COHORT {

    tag "cohort"

    input:
        val token
        tuple path(manifest), val(use_demo)

    output:
        val true

    script:
    """
    cd "${projectDir}"

    if [ "${use_demo}" = "true" ]; then
        bash bin/cohort.sh ${manifest} --demo
    else
        bash bin/cohort.sh ${manifest}
    fi
    """
}

process COHORT_FILTER {

    tag "cohort-filter"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    python bin/cohortFilter.py
    """
}

process SNP_MATRIX {

    tag "snp-matrix"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    python bin/snpMatrix.py
    """
}

process TRANSMISSION {

    tag "transmission"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    python bin/transmission.py
    """
}

process IQTREE {

    tag "iqtree"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    bash bin/iqtree.sh
    """
}

/*
 * ============================================================
 * BLOCO 3 — PER BIOSAMPLE
 * Inicia apenas após sincronização global do BLOCO 2
 * ============================================================
 */

process MIX_INFECTION {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    python bin/mixInfection.py ${biosample}
    """
}


process RESISTANCE_TARGET {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    python bin/resistanceTarget.py ${biosample}
    """
}


process RESISTANCE_REPORT {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    python bin/resistanceReport.py ${biosample}
    """
}

/*
 * ============================================================
 * BLOCO  — COHORT e PER BIOSAMPLE
 * Inicia apenas após sincronização global do BLOCO 3
 * ============================================================
 */
 
process RESISTANCE_SUMMARY {

    tag "resistance-summary"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    python bin/resistanceSummary.py
    """
}

process QC_SUMMARY {

    tag "qc-summary"

    input:
        val token

    output:
        val true

    script:
    """
    cd "${projectDir}"
    python bin/qcSummary.py
    """
}

process CLINICAL_REPORT {

    tag { biosample }

    input:
        val biosample

    output:
        val biosample

    script:
    """
    cd "${projectDir}"
    python bin/clinicalReport.py ${biosample}
    """
}

workflow {

    /*
     * ============================================================
     * INIT (SEQUENTIAL) — SEMPRE EXECUTA
     * ============================================================
     */

    init_done = Channel.value(true)
        | KAIJU_DB
        | OMS_CATALOG
        | BWA_REF
        | GATK_DICT
        | SNPEFF_DB

    /*
     * ============================================================
     * MANIFEST — SEMPRE EXECUTA
     * ============================================================
     */

    manifest_ch = MAKE_MANIFEST_VALIDATE(
        init_done,
        file(params.input_table),
        file(params.reads_dir)
    )

    /*
     * ============================================================
     * SAMPLES CHANNEL
     * ============================================================
     */

    samples_ch = manifest_ch
        .splitCsv(header:true, sep:'\t')
        .map { row -> row.biosample }


    /*
     * ============================================================
     * EXCLUSÕES PERMITIDAS (SOMENTE FULL PIPELINE)
     * ============================================================
     */

    def excluded = params.exclude ?
        params.exclude.split(',')*.trim() :
        []

    def allowed_exclusions = ['kaiju','iqtree','transmission','clinical_report']

    excluded.each {
        if (!allowed_exclusions.contains(it)) {
            error "Exclusão inválida: ${it}. Permitidos: ${allowed_exclusions}"
        }
    }


    /*
     * ============================================================
     * MODO MÓDULO (EXECUÇÃO ISOLADA)
     * ============================================================
     */

    if (params.module) {

        switch (params.module) {

            case 'fastqc': FASTQC(samples_ch); break
            case 'trimmomatic': TRIMMOMATIC(samples_ch); break
            case 'kaiju': KAIJU(samples_ch); break
            case 'bwa': BWA(samples_ch); break
            case 'delly': DELLY(samples_ch); break
            case 'lofreq': LOFREQ(samples_ch); break
            case 'gatk_gvcf': GATK_GVCF(samples_ch); break
            case 'gatk_vcf': GATK_VCF(samples_ch); break
            case 'norm': NORM(samples_ch); break
            case 'tbdr_rcov': TBDR_RCOV(samples_ch); break
            case 'lineage': LINEAGE(samples_ch); break
            case 'ntm_filter': NTM_FILTER(samples_ch); break
            case 'snpeff': SNPEFF(samples_ch); break

            case 'cohort':
                cohort_input_ch = manifest_ch
                    .map { file -> tuple(file, params.demo) }
                COHORT(Channel.value(true), cohort_input_ch)
                break

            case 'cohort_filter': COHORT_FILTER(Channel.value(true)); break
            case 'snp_matrix': SNP_MATRIX(Channel.value(true)); break
            case 'transmission': TRANSMISSION(Channel.value(true)); break
            case 'iqtree': IQTREE(Channel.value(true)); break
            case 'mix_infection': MIX_INFECTION(samples_ch); break
            case 'resistance_target': RESISTANCE_TARGET(samples_ch); break
            case 'resistance_report': RESISTANCE_REPORT(samples_ch); break
            case 'resistance_summary': RESISTANCE_SUMMARY(Channel.value(true)); break
            case 'qc_summary': QC_SUMMARY(Channel.value(true)); break
            case 'clinical_report': CLINICAL_REPORT(samples_ch); break

            default:
                error "Módulo inválido: ${params.module}"
        }

    } else {

        /*
         * ============================================================
         * PIPELINE COMPLETA
         * ============================================================
         */

        // BLOCO 1

        fastqc_ch = FASTQC(samples_ch)
        trimmomatic_ch = TRIMMOMATIC(samples_ch)

        if (!excluded.contains('kaiju')) {
            kaiju_ch = KAIJU(trimmomatic_ch)
        }

        bwa_ch   = BWA(trimmomatic_ch)

        delly_ch      = DELLY(bwa_ch)
        lofreq_ch     = LOFREQ(bwa_ch)
        gatk_gvcf_ch  = GATK_GVCF(bwa_ch)
        gatk_vcf_ch   = GATK_VCF(bwa_ch)

        norm_ch = NORM(gatk_vcf_ch)
        tbdr_rcov_ch = TBDR_RCOV(gatk_gvcf_ch)

        lineage_ch = LINEAGE(norm_ch)

        ntm_input_ch = lofreq_ch.join(gatk_gvcf_ch)
        ntm_filter_ch = NTM_FILTER(ntm_input_ch)

        snpeff_input_ch = lofreq_ch
            .join(gatk_gvcf_ch)
            .join(gatk_vcf_ch)
            .join(norm_ch)

        snpeff_ch = SNPEFF(snpeff_input_ch)
        
          /*
         * ============================================================
         * BLOCO 1 — SINCRONIZAÇÃO CORRETA
         * ============================================================
         */

        bloco1_barrier_channels = []

        bloco1_barrier_channels << snpeff_ch

        if (!excluded.contains('kaiju')) {
            bloco1_barrier_channels << kaiju_ch
        }

        if (bloco1_barrier_channels.size() == 2) {

            bloco1_sync = bloco1_barrier_channels[0]
                .join(bloco1_barrier_channels[1])
                .collect()
                .map { true }

        } else {

            bloco1_sync = bloco1_barrier_channels[0]
                .collect()
                .map { true }
        }


        
        // BLOCO 2

        cohort_input_ch = manifest_ch
            .map { file -> tuple(file, params.demo) }

        cohort_ch = COHORT(bloco1_sync, cohort_input_ch)

        cohort_filter_ch = COHORT_FILTER(cohort_ch)
        snp_matrix_ch    = SNP_MATRIX(cohort_filter_ch)

        def barrier_channels = []

        if (!excluded.contains('transmission')) {
            transmission_ch = TRANSMISSION(snp_matrix_ch)
            barrier_channels << transmission_ch
        }

        if (!excluded.contains('iqtree')) {
            iqtree_ch = IQTREE(snp_matrix_ch)
            barrier_channels << iqtree_ch
        }

        if (barrier_channels.size() == 2) {

            bloco2_sync = barrier_channels[0]
                .join(barrier_channels[1])
                .collect()
                .map { true }

        } else if (barrier_channels.size() == 1) {

            bloco2_sync = barrier_channels[0]
                .collect()
                .map { true }

        } else {

            bloco2_sync = snp_matrix_ch
                .collect()
                .map { true }
        }


        // BLOCO 3

        bloco3_samples_ch = samples_ch
            .combine(bloco2_sync)
            .map { biosample, _ -> biosample }

        mix_infection_ch      = MIX_INFECTION(bloco3_samples_ch)
        resistance_target_ch  = RESISTANCE_TARGET(bloco3_samples_ch)
        resistance_report_ch  = RESISTANCE_REPORT(resistance_target_ch)

        bloco3_sync = resistance_report_ch.collect().map { true }


        // BLOCO 4

        resistance_summary_ch = RESISTANCE_SUMMARY(bloco3_sync)
        qc_summary_ch         = QC_SUMMARY(bloco3_sync)

        bloco4_sync = resistance_summary_ch
            .join(qc_summary_ch)
            .collect()
            .map { true }

        clinical_samples_ch = samples_ch
            .combine(bloco4_sync)
            .map { biosample, _ -> biosample }

        if (!excluded.contains('clinical_report')) {
            clinical_report_ch = CLINICAL_REPORT(clinical_samples_ch)
        }
    }
}


