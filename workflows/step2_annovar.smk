chr_name_list = [c+1 for c in range(22)]
pop_name_list = ['ACB', 'BEB', 'CEU', 'CHS', 'ESN', 'GBR', 'GWD', 'ITU', 'KHV', 'MSL', 'PEL', 'PUR', 'TSI', 'ASW', 'CDX', 'CHB', 'CLM', 'FIN', 'GIH', 'IBS', 'JPT', 'LWK', 'MXL', 'PJL', 'STU', 'YRI']

rule all:
    input:
        expand("results/exonic_data/ALL.{mut_type}.vcf.gz", mut_type=['synonymous', 'nonsynonymous']),
        expand("results/populations/{pop_name}.list", pop_name=pop_name_list),


rule extract_pop_info:
    input:
        ped = "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
    output:
        expand("results/populations/{pop_name}.list", pop_name=pop_name_list)
    shell:
        """
        grep -w 'ACB' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/ACB.list
        grep -w 'BEB' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/BEB.list
        grep -w 'CEU' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/CEU.list
        grep -w 'CHS' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/CHS.list
        grep -w 'ESN' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/ESN.list
        grep -w 'GBR' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/GBR.list
        grep -w 'GWD' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/GWD.list
        grep -w 'ITU' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/ITU.list
        grep -w 'KHV' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/KHV.list
        grep -w 'MSL' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/MSL.list
        grep -w 'PEL' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/PEL.list
        grep -w 'PUR' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/PUR.list
        grep -w 'TSI' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/TSI.list
        grep -w 'ASW' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/ASW.list
        grep -w 'CDX' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/CDX.list
        grep -w 'CHB' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/CHB.list
        grep -w 'CLM' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/CLM.list
        grep -w 'FIN' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/FIN.list
        grep -w 'GIH' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/GIH.list
        grep -w 'IBS' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/IBS.list
        grep -w 'JPT' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/JPT.list
        grep -w 'LWK' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/LWK.list
        grep -w 'MXL' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/MXL.list
        grep -w 'PJL' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/PJL.list
        grep -w 'STU' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/STU.list
        grep -w 'YRI' {input.ped} | awk '{{print $1"\\t"$2}}' > results/populations/YRI.list
        """


rule extract_biallelic_snps:
    input:
        vcf = "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr_name}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    output:
        vcf = "results/genotypes/ALL.chr{chr_name}.phase3_shapeit2_mvncall_integrated_v5b.20130502.biallelic.snps.genotypes.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule download_annovar_db:
    input:
    output:
        avsnp150 = "ext/annovar/humandb/hg19_avsnp150.txt",
        dbnsfp42c = "ext/annovar/humandb/hg19_dbnsfp42c.txt",
    shell:
        """
        ext/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar avsnp150 ext/annovar/humandb/
        ext/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar dbnsfp42c ext/annovar/humandb/
        """


rule run_annovar:
    input:
        vcf = rules.extract_biallelic_snps.output.vcf,
        avsnp150 = rules.download_annovar_db.output.avsnp150,
        dbnsfp42c = rules.download_annovar_db.output.dbnsfp42c,
    output:
        vcf = "results/genotypes/ALL.chr{chr_name}.annotated.biallelic.snps.hg19_multianno.vcf",
        txt = "results/genotypes/ALL.chr{chr_name}.annotated.biallelic.snps.hg19_multianno.txt",
    resources:
        time=2880, cpus=8, mem_gb=64,
    params:
        output_prefix = "results/genotypes/ALL.chr{chr_name}.annotated.biallelic.snps",
    shell:
        """
        ext/annovar/table_annovar.pl {input.vcf} ext/annovar/humandb/ -buildver hg19 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . -vcfinput --thread {resources.cpus}
        """


rule compress:
    input:
        vcf = rules.run_annovar.output.vcf,
    output:
        vcf = "results/genotypes/ALL.chr{chr_name}.annotated.biallelic.snps.hg19_multianno.vcf.gz",
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {input.vcf}
        """


rule extract_exonic_snps:
    input:
        vcf = rules.extract_biallelic_snps.output.vcf,
        file = rules.run_annovar.output.txt,
    output:
        synonymous_vcf = "results/exonic_data/ALL.chr{chr_name}.synonymous.vcf.gz",
        nonsynonymous_vcf = "results/exonic_data/ALL.chr{chr_name}.nonsynonymous.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -R <(grep -w synonymous {input.file} | awk '{{print $1"\\t"$2}}') | bgzip -c > {output.synonymous_vcf}
        bcftools view {input.vcf} -R <(grep -w nonsynonymous {input.file} | awk '{{print $1"\\t"$2}}') | bgzip -c > {output.nonsynonymous_vcf}
        tabix -p vcf {output.synonymous_vcf}
        tabix -p vcf {output.nonsynonymous_vcf}
        """


rule concat_files:
    input:
        synonymous_vcfs = expand("results/exonic_data/ALL.chr{chr_name}.synonymous.vcf.gz", chr_name=chr_name_list),
        nonsynonymous_vcfs = expand("results/exonic_data/ALL.chr{chr_name}.nonsynonymous.vcf.gz", chr_name=chr_name_list),
    output:
        synonymous_vcf = "results/exonic_data/ALL.synonymous.vcf.gz",
        nonsynonymous_vcf = "results/exonic_data/ALL.nonsynonymous.vcf.gz",
    shell:
        """
        bcftools concat {input.synonymous_vcfs} | bgzip -c > {output.synonymous_vcf}
        bcftools concat {input.nonsynonymous_vcfs} | bgzip -c > {output.nonsynonymous_vcf}
        tabix -p vcf {output.synonymous_vcf}
        tabix -p vcf {output.nonsynonymous_vcf}
        """
