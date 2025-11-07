awk -F'\t' '$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, a);
    match($9, /description "([^"]+)"/, b);
    gene_id = (a[1] != "") ? a[1] : "NA";
    description = (b[1] != "") ? b[1] : "NA";
    print $1"\t"gene_id"\t"description"\t"$4"\t"$5"\t"$7
}' GCF_019914945.1_TyphaL0001v2_genomic.gtf > gene_info_with_description.tsv

