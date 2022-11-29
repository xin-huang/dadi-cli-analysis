rule all:
    input:
        ".downloaded",
        

rule download:
    input:
    output:
        ".downloaded",
    shell:
        """
        wget -rc -np -l 1 -R "index.html*" -e robots=off https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
        touch .downloaded
        """
