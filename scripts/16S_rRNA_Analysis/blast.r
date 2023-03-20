#Running blast on the server
blastn = "/opt/apps/blast/2.10.1+/bin/blastn"
blast_db = "./blast_dir/16S_ribosomal_RNA"
input = "./unassigned.fasta"
evalue = 9.6e-6
format = 6
max_target = 1
colnames <- c("qseqid",
              "sseqid",
              "evalue",
              "bitscore",
              "sgi",
              "sacc")
blast_out <- system2("/opt/apps/blast/2.10.1+/bin/blastn",
                     args = c("-db", blast_db,
                              "-query", input,
                              "-outfmt", format,
                              "-evalue", evalue,
                              "-max_target_seqs", max_target,
                              "-ungapped"),
                     wait = TRUE,
                     stdout = TRUE) %>%
  as_tibble() %>%
  separate(col = value,
           into = colnames,
           sep = "\t",
           convert = TRUE)

