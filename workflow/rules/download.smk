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
        "results/data/1KG/.downloaded",
    params:
        dir = "results/data/1KG"
    log:
        "logs/download/download_1KG.log"
    shell:
        """
        cd {params.dir}
        wget -rc -np -l 1 -R "index.html*" -e robots=off https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
        touch .downloaded
        """
