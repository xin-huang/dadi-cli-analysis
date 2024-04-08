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


rule download_1KG:
    input:
    output:
        "results/1KG/data/.downloaded",
    params:
        dir = "results/1KG/data",
        dl_dir = "ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    log:
        "logs/download/download.1KG.log"
    shell:
        """
        cd {params.dir}
        wget -rc -np -l 1 -R "index.html*" -e robots=off https://{params.dl_dir}
        mv {params.dl_dir}* .
        rm -r {params.dl_dir}
        touch .downloaded
        """


rule download_annovar_db:
    input:
    output:
        "resources/annovar/humandb/{ref_gene}_{database}.txt",
    log:
        "logs/download/download.annovar.db.{ref_gene}.{database}.log"
    shell:
        """
        resources/annovar/annotate_variation.pl -downdb -buildver {wildcards.refgen} \
             -webfrom annovar {wildcards.database} resources/annovar/humandb/
        """
