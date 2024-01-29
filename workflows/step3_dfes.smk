population_list = [
    'ACB',
    'ASW',
    'BEB',
    'CDX',
    'CEU',
    'CHB',
    'CHS',
    'CLM',
    'ESN',
    'FIN',
    'GBR',
    'GIH',
    'GWD',
    'IBS',
    'ITU',
    'JPT',
    'KHV',
    'LWK',
    'MSL',
    'MXL',
    'PEL',
    'PJL',
    'PUR',
    'STU',
    'TSI',
    'YRI',
]

sample_size_list = {
    'ACB': 96*2,
    'ASW': 61*2,
    'BEB': 86*2,
    'CDX': 93*2,
    'CEU': 99*2,
    'CHB': 103*2,
    'CHS': 105*2,
    'CLM': 94*2,
    'ESN': 99*2,
    'FIN': 99*2,
    'GBR': 91*2,
    'GIH': 103*2,
    'GWD': 113*2,
    'IBS': 107*2,
    'ITU': 102*2,
    'JPT': 104*2,
    'KHV': 99*2,
    'LWK': 99*2,
    'MSL': 85*2,
    'MXL': 64*2,
    'PEL': 85*2,
    'PJL': 96*2,
    'PUR': 104*2,
    'STU': 102*2,
    'TSI': 107*2,
    'YRI': 108*2,
}

model_list = ['two_epoch']
model_params_list = {
    'two_epoch': {
        'p0': '5 0.5 0.5',
        'ubounds': '100 1 1',
        'lbounds': '10e-3 10e-3 0',
    },
}
dfe_list = ['gamma', 'lognormal']
dfe_params_list = {
    'gamma': {
        'p0': '50 25000 0.5',
        'ubounds': '100 50000 1',
        'lbounds': '10e-3 10e-3 0',
    },
    'lognormal': {
        'p0': '5 5 0.5',
        'ubounds': '100 100 1',
        'lbounds': '10e-3 10e-3 0',
    },
}


rule all:
    input:
        expand("results/dfes/{population}/unfolded/StatDFE/{population}.{model}.{dfe}.godambe.ci", population=population_list, model=model_list, dfe=dfe_list),
        "results/dfes/all.res.txt",


rule extract_samples:
    input:
        synonymous_vcf = "results/exonic_data/ALL.synonymous.vcf.gz",
        nonsynonymous_vcf = "results/exonic_data/ALL.nonsynonymous.vcf.gz",
        samples = "results/populations/{population}.list",
    output:
        synonymous_vcf = "results/dfes/{population}/inputs/{population}.synonymous.biallelic.snps.vcf.gz",
        nonsynonymous_vcf = "results/dfes/{population}/inputs/{population}.nonsynonymous.biallelic.snps.vcf.gz",
        info = "results/dfes/{population}/inputs/{population}.info",
    shell:
        """
        bcftools view {input.synonymous_vcf} -S <(awk '{{print $1}}' {input.samples}) --force-samples -v snps -m 2 -M 2 | bgzip -c > {output.synonymous_vcf}
        bcftools view {input.nonsynonymous_vcf} -S <(awk '{{print $1}}' {input.samples}) --force-samples -v snps -m 2 -M 2 | bgzip -c > {output.nonsynonymous_vcf}
        tabix -p vcf {output.synonymous_vcf}
        tabix -p vcf {output.nonsynonymous_vcf}
        bcftools query -l {output.synonymous_vcf} | awk -v pop={wildcards.population} '{{print $1"\\t"pop}}' > {output.info}
        """


rule add_AA:
    input:
        synonymous_vcf = rules.extract_samples.output.synonymous_vcf,
        nonsynonymous_vcf = rules.extract_samples.output.nonsynonymous_vcf,
    output:
        synonymous_AA_vcf = "results/dfes/{population}/inputs/{population}.synonymous.biallelic.snps.AA.vcf.gz",
        nonsynonymous_AA_vcf = "results/dfes/{population}/inputs/{population}.nonsynonymous.biallelic.snps.AA.vcf.gz",
    params:
        synonymous_AA_bed = "results/dfes/{population}/inputs/{population}.synonymous.AA.bed.gz",
        nonsynonymous_AA_bed = "results/dfes/{population}/inputs/{population}.nonsynonymous.AA.bed.gz",
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%INFO/AA\\n" {input.synonymous_vcf} | sed 's/|//g' | sed 's/a/A/' | sed 's/c/C/' | sed 's/t/T/' | sed 's/g/G/' | awk '{{print $1"\\t"$2-1"\\t"$2"\\t"$3}}' | bgzip -c > {params.synonymous_AA_bed}
        bcftools query -f "%CHROM\\t%POS\\t%INFO/AA\\n" {input.nonsynonymous_vcf} | sed 's/|//g' | sed 's/a/A/' | sed 's/c/C/' | sed 's/t/T/' | sed 's/g/G/' | awk '{{print $1"\\t"$2-1"\\t"$2"\\t"$3}}' | bgzip -c > {params.nonsynonymous_AA_bed}
        tabix -p bed {params.synonymous_AA_bed}
        tabix -p bed {params.nonsynonymous_AA_bed}
        bcftools annotate -a {params.synonymous_AA_bed} -h <(echo "##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">") -c CHROM,FROM,TO,AA {input.synonymous_vcf} | bgzip -c > {output.synonymous_AA_vcf}
        bcftools annotate -a {params.nonsynonymous_AA_bed} -h <(echo "##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">") -c CHROM,FROM,TO,AA {input.nonsynonymous_vcf} | bgzip -c > {output.nonsynonymous_AA_vcf}
        """


rule generate_fs:
    input:
        synonymous_vcf = rules.add_AA.output.synonymous_AA_vcf,
        nonsynonymous_vcf = rules.add_AA.output.nonsynonymous_AA_vcf,
        info = rules.extract_samples.output.info,
    output:
        synonymous_fs = "results/dfes/{population}/unfolded/{population}.synonymous.unfolded.fs",
        nonsynonymous_fs = "results/dfes/{population}/unfolded/{population}.nonsynonymous.unfolded.fs",
    params:
        sample_size = lambda wildcards: sample_size_list[wildcards.population]
    shell:
        """
        dadi-cli GenerateFs --vcf {input.synonymous_vcf} --pop-ids {wildcards.population} --pop-info {input.info} --projections {params.sample_size} --output {output.synonymous_fs} --polarized
        dadi-cli GenerateFs --vcf {input.nonsynonymous_vcf} --pop-ids {wildcards.population} --pop-info {input.info} --projections {params.sample_size} --output {output.nonsynonymous_fs} --polarized
        """


rule infer_dm:
    input:
        fs = rules.generate_fs.output.synonymous_fs,
    output:
        bestfit = "results/dfes/{population}/unfolded/InferDM/{population}.{model}.InferDM.bestfits",
    resources:
        cpus = 8,
    params:
        model = "{model}",
        p0 = lambda wildcards: model_params_list[wildcards.model]['p0'],
        ubounds = lambda wildcards: model_params_list[wildcards.model]['ubounds'],
        lbounds = lambda wildcards: model_params_list[wildcards.model]['lbounds'],
        grid_size = '300 400 500',
        output_prefix = "results/dfes/{population}/unfolded/InferDM/{population}.{model}",
        optimizations = 100,
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {params.model} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --output-prefix {params.output_prefix} --grids {params.grid_size} --cpus {resources.cpus} --force-convergence {params.optimizations}
        """


rule generate_cache:
    input:
        bestfit = rules.infer_dm.output.bestfit,
    output:
        cache = "results/dfes/{population}/unfolded/InferDFE/{population}.{model}.spectra.bpkl",
    resources:
        cpus = 8,
    params:
        model = "{model}_sel",
        sample_size = lambda wildcards: sample_size_list[wildcards.population],
        grid_size = "800 1000 1200",
        gamma_pts = 2000,
    shell:
        """
        dadi-cli GenerateCache --model {params.model} --demo-popt {input.bestfit} --sample-size {params.sample_size} --output {output.cache} --cpus {resources.cpus} --grids {params.grid_size} --gamma-pts {params.gamma_pts}
        """


rule infer_dfe:
    input:
        fs = rules.generate_fs.output.nonsynonymous_fs,
        cache = rules.generate_cache.output.cache,
        dm_bestfit = rules.infer_dm.output.bestfit,
    output:
        bestfit = "results/dfes/{population}/unfolded/InferDFE/{population}.{model}.{dfe}.InferDFE.bestfits",
    resources:
        cpus = 8,
    params:
        dfe = "{dfe}",
        p0 = lambda wildcards: dfe_params_list[wildcards.dfe]['p0'],
        ubounds = lambda wildcards: dfe_params_list[wildcards.dfe]['ubounds'],
        lbounds = lambda wildcards: dfe_params_list[wildcards.dfe]['lbounds'],
        ratio = 2.31,
        output_prefix = "results/dfes/{population}/unfolded/InferDFE/{population}.{model}.{dfe}", 
        optimizations = 100,
    shell:
        """
        dadi-cli InferDFE --fs {input.fs} --cache1d {input.cache} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf1d {params.dfe} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --ratio {params.ratio} --cpus {resources.cpus} --force-convergence {params.optimizations}
        """


rule godambe_ci:
    input:
        synonymous_vcf = rules.add_AA.output.synonymous_AA_vcf,
        nonsynonymous_vcf = rules.add_AA.output.nonsynonymous_AA_vcf,
        popinfo = rules.extract_samples.output.info,
        nonsynonymous_fs = rules.generate_fs.output.nonsynonymous_fs,
        cache = rules.generate_cache.output.cache,
        dfe_bestfit = rules.infer_dfe.output.bestfit,
    output:
        dfe_godambe_ci = "results/dfes/{population}/unfolded/StatDFE/{population}.{model}.{dfe}.godambe.ci",
    params:
        synonymous_dir = "results/dfes/{population}/unfolded/{population}_bootstrapping_syn",
        nonsynonymous_dir = "results/dfes/{population}/unfolded/{population}_bootstrapping_non",
        synonymous_output_prefix = "results/dfes/{population}/unfolded/{population}_bootstrapping_syn/{population}.synonymous.unfolded",
        nonsynonymous_output_prefix = "results/dfes/{population}/unfolded/{population}_bootstrapping_non/{population}.nonsynonymous.unfolded",
        sample_size = lambda wildcards: sample_size_list[wildcards.population],
    shell:
        """
        [ -d {params.synonymous_dir} ] || mkdir {params.synonymous_dir}
        [ -d {params.nonsynonymous_dir} ] || mkdir {params.nonsynonymous_dir}
        dadi-cli GenerateFs --vcf {input.synonymous_vcf} --pop-info {input.popinfo} --pop-ids {wildcards.population} --projections {params.sample_size} --polarized --bootstrap 100 --chunk-size 1000000 --output {params.synonymous_output_prefix}
        dadi-cli GenerateFs --vcf {input.nonsynonymous_vcf} --pop-info {input.popinfo} --pop-ids {wildcards.population} --projections {params.sample_size} --polarized --bootstrap 100 --chunk-size 1000000 --output {params.nonsynonymous_output_prefix}
        dadi-cli StatDFE --fs {input.nonsynonymous_fs} --dfe-popt {input.dfe_bestfit} --cache1d {input.cache} --pdf1d {wildcards.dfe} --bootstrapping-nonsynonymous-dir {params.synonymous_dir} --bootstrapping-synonymous-dir {params.nonsynonymous_dir} --output {output.dfe_godambe_ci}
        """


rule single_pop_summary:
    input:
        dfe_bestfit = rules.infer_dfe.output.bestfit,
        dfe_godambe_ci = rules.godambe_ci.output.dfe_godambe_ci,
    output:
        res = "results/dfes/{population}/unfolded/{population}.{model}.{dfe}.res.txt",
    shell:
        """
        paste <(grep Converged {input.dfe_bestfit} -A 2 | tail -1) <(grep "step size 0.001" {input.dfe_godambe_ci} -A 1 | tail -1 | sed 's/\\[//' | sed 's/\\]//' | awk '{{print $(NF-2)"\\t"$(NF-1)"\\t"$NF}}') <(grep "step size 0.001" {input.dfe_godambe_ci} -A 2 | tail -1 | sed 's/\\[//' | sed 's/\\]//' | awk '{{print $(NF-2)"\\t"$(NF-1)"\\t"$NF}}') | awk -v pop={wildcards.population} -v demog={wildcards.model} -v dfe={wildcards.dfe} '{{print pop"\\t"demog"\\t"dfe"\\t"$0}}' > {output.res}
        """


rule all_pops_summary:
    input:
        res = expand("results/dfes/{population}/unfolded/{population}.{model}.{dfe}.res.txt", population=population_list, model=model_list, dfe=dfe_list),
    output:
        res = "results/dfes/all.res.txt",
    shell:
        """
        cat {input.res} | sed '1iPop\\tDemog\\tDFE\\tLog(likelihood)\\tDFE_param1\\tDFE_param2\\tMisid\\tTheta\\tDFE_param1_lb\\tDFE_param2_lb\\tMisid_lb\\tDFE_param1_ub\\tDFE_param2_ub\\tMisid_ub' > {output.res}
        """

