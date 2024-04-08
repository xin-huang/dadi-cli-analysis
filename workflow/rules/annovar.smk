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
        vcf = get_input_vcf,
    output:
        vcf = "results/{dataset}/biallelic_snps/{output_prefix}.chr{chr_name}.vcf.gz",
    log:
        "logs/extract_biallelic_snps/{dataset}/{output_prefix}.{chr_name}.log",
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule run_annovar:
    input:
        vcf = rules.extract_biallelic_snps.output.vcf,
        avsnp150 = "resources/annovar/humandb/{ref_gene}_avsnp150.txt",
        dbnsfp42c = "resources/annovar/humandb/{ref_gene}_dbnsfp42c.txt",
    output:
        vcf = "results/{dataset}/annotated_snps/{output_prefix}.chr{chr_name}.annotated.biallelic.snps.{ref_gene}_multianno.vcf.gz",
        txt = "results/{dataset}/annotated_snps/{output_prefix}.chr{chr_name}.annotated.biallelic.snps.{ref_gene}_multianno.txt",
    log:
        "logs/run_annovar/{dataset}/{output_prefix}.chr{chr_name}.{ref_gene}.log",
    resources:
        time = 2880, 
        cpus = 8, 
        mem_gb = 64,
    params:
        vcf = "results/{dataset}/annotated_snps/{output_prefix}.chr{chr_name}.annotated.biallelic.snps.{ref_gene}_multianno.vcf",
        output_prefix = "results/{dataset}/annotated_snps/{output_prefix}.chr{chr_name}.annotated.biallelic.snps",
    shell:
        """
        resources/annovar/table_annovar.pl {input.vcf} resources/annovar/humandb/ -buildver {wildcards.ref_gene} -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . -vcfinput --thread {resources.cpus}
        bgzip -c {params.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {params.vcf}
        """


rule extract_exonic_snps:
    input:
        vcf = rules.run_annovar.output.vcf,
        file = rules.run_annovar.output.txt,
    output:
        vcf = "results/{dataset}/exonic_data/{output_prefix}.chr{chr_name}.{ref_gene}.{mut_type}.vcf.gz",
    log:
        "logs/extract_exonic_snps/{dataset}/{output_prefix}.{chr_name}.{ref_gene}.{mut_type}.log"
    shell:
        """
        bcftools view {input.vcf} -R <(grep -w {wildcards.mut_type} {input.file} | awk '{{print $1"\\t"$2}}') | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule concat_files:
    input:
        vcfs = expand("results/{dataset}/exonic_data/{output_prefix}.chr{chr_name}.{ref_gene}.{mut_type}.vcf.gz", 
                      chr_name=chr_name_list, allow_missing=True),
    output:
        vcf = "results/{dataset}/exonic_data/{output_prefix}.{ref_gene}.{mut_type}.vcf.gz",
    log:
        "logs/concat_files/{dataset}/{output_prefix}.{ref_gene}.{mut_type}.log"
    shell:
        """
        bcftools concat {input.vcfs} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
