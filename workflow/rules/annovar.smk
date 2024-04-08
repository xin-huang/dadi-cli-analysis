# GNU General Public License v3.0
# Copyright 2024 Xin Huang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see 
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


rule extract_biallelic_snps:
    input:
        vcf = "results/{dataset}/data/ALL.chr{chr_name}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    output:
        vcf = "results/{dataset}/biallelic_snps/ALL.chr{chr_name}.phase3_shapeit2_mvncall_integrated_v5b.20130502.biallelic.snps.genotypes.vcf.gz",
    log:
        "logs/extract_biallelic_snps/{dataset}/{chr_name}.log",
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule run_annovar:
    input:
        vcf = rules.extract_biallelic_snps.output.vcf,
        avsnp150 = rules.download_annovar_db.output.avsnp150,
        dbnsfp42c = rules.download_annovar_db.output.dbnsfp42c,
    output:
        vcf = "results/{dataset}/annotated_snps/ALL.chr{chr_name}.annotated.biallelic.snps.hg19_multianno.vcf.gz",
        txt = "results/{dataset}/annotated_snps/ALL.chr{chr_name}.annotated.biallelic.snps.hg19_multianno.txt",
    log:
        "logs/run_annovar/{dataset}/{chr_name}.log",
    resources:
        time = 2880, 
        cpus = 8, 
        mem_gb = 64,
    params:
        vcf = "results/{dataset}/annotated_snps/ALL.chr{chr_name}.annotated.biallelic.snps.hg19_multianno.vcf",
        output_prefix = "results/{dataset}/annotated_snps/ALL.chr{chr_name}.annotated.biallelic.snps",
    shell:
        """
        resources/annovar/table_annovar.pl {input.vcf} resources/annovar/humandb/ \
            -buildver hg19 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c \
            -operation g,f,f -nastring . -vcfinput --thread {resources.cpus}
        bgzip -c {params.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {params.vcf}
        """


rule extract_exonic_snps:
    input:
        vcf = rules.compress.output.vcf,
        file = rules.run_annovar.output.txt,
    output:
        synonymous_vcf = "results/{dataset}/exonic_data/ALL.chr{chr_name}.synonymous.vcf.gz",
        nonsynonymous_vcf = "results/{dataset}/exonic_data/ALL.chr{chr_name}.nonsynonymous.vcf.gz",
    log:
        "logs/extract_exonic_snps/{dataset}/{chr_name}.log"
    shell:
        """
        bcftools view {input.vcf} -R <(grep -w synonymous {input.file} | awk '{{print $1"\\t"$2}}') | bgzip -c > {output.synonymous_vcf}
        bcftools view {input.vcf} -R <(grep -w nonsynonymous {input.file} | awk '{{print $1"\\t"$2}}') | bgzip -c > {output.nonsynonymous_vcf}
        tabix -p vcf {output.synonymous_vcf}
        tabix -p vcf {output.nonsynonymous_vcf}
        """


rule concat_files:
    input:
        synonymous_vcfs = expand("results/{dataset}/exonic_data/ALL.chr{chr_name}.synonymous.vcf.gz", chr_name=chr_name_list),
        nonsynonymous_vcfs = expand("results/{dataset}/exonic_data/ALL.chr{chr_name}.nonsynonymous.vcf.gz", chr_name=chr_name_list),
    output:
        synonymous_vcf = "results/{dataset}/exonic_data/ALL.synonymous.vcf.gz",
        nonsynonymous_vcf = "results/{dataset}/exonic_data/ALL.nonsynonymous.vcf.gz",
    log:
        "logs/concat_files/{dataset}/concat.log"
    shell:
        """
        bcftools concat {input.synonymous_vcfs} | bgzip -c > {output.synonymous_vcf}
        bcftools concat {input.nonsynonymous_vcfs} | bgzip -c > {output.nonsynonymous_vcf}
        tabix -p vcf {output.synonymous_vcf}
        tabix -p vcf {output.nonsynonymous_vcf}
        """
