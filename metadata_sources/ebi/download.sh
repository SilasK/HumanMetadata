metadata_file='genomes-all_metadata.tsv'

base_uhgg="http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0"

wget ${metadata_file}  ${base_uhgg}/${metadata_file}
