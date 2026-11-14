# ============================================================
# MOLECULAR MECHANISMS OF MUSCLE CONTRACTILITY
# Gene: CSRP3
# Assignment No.: 20
# Variant: NM_003476.5:c.172T>G (p.Cys58Gly)
# ClinVar Variation ID: 8777
# ============================================================


# STEP 1. SETUP
# Install packages only if they are not yet installed.

if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(rentrez)
library(ggplot2)

# Assignment information
gene_symbol <- "CSRP3"
transcript_accession <- "NM_003476.5"
clinvar_id <- "8777"
variant_protein <- "p.Cys58Gly"

# Create output folders
for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}


# STEP 2. RETRIEVE HUMAN NCBI GENE RECORD

hits <- entrez_search(
  db = "gene",
  term = paste0(
    gene_symbol,
    "[Gene Name] AND txid9606[Organism]"
  ),
  retmax = 20
)

if (length(hits$ids) == 0) {
  stop("Gene not found")
}

# Verify the exact gene symbol
gene_id <- NA_character_

for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id
    gene_info <- rec
    break
  }
}

if (is.na(gene_id)) {
  stop("Exact human gene not found")
}

print(gene_info)

# Save gene information
writeLines(
  c(
    paste("Gene symbol:", gene_info$name),
    paste("NCBI Gene ID:", gene_id),
    paste("Official name:", gene_info$description),
    paste("Chromosome:", gene_info$chromosome)
  ),
  "results/gene_summary.txt"
)


# STEP 3. RETRIEVE MATCHING PROTEIN SEQUENCE
# Retrieve the exact assigned NM_ transcript.

gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

# Extract the RefSeq NP_ protein accession
pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'

fragment <- regmatches(
  gb,
  regexpr(pat, gb)
)

if (length(fragment) == 0) {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub(
  '^/protein_id="|"$',
  "",
  fragment
)

# Retrieve the protein FASTA sequence
protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Invalid FASTA")
}

# Save FASTA file
fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)

cat("Matched protein:", protein_accession, "\n")

# Save protein accession
writeLines(
  paste("NP accession:", protein_accession),
  "results/protein_accession.txt"
)


# STEP 4. ANALYZE PROTEIN LENGTH AND AMINO-ACID COMPOSITION

# Read the saved FASTA
fasta_lines <- readLines(fasta_file)

sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

aa <- strsplit(toupper(sequence), "")[[1]]

# Calculate protein length
protein_length <- length(aa)

# Count the 20 standard amino acids
standard <- strsplit(
  "ACDEFGHIKLMNPQRSTVWY",
  ""
)[[1]]

counts <- table(
  factor(aa, levels = standard)
)

aa_results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

# Save amino-acid composition table
write.csv(
  aa_results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

# Display protein length and top three amino acids
cat(
  "Protein length:",
  protein_length,
  "amino acids\n"
)

print(
  head(
    aa_results[order(-aa_results$Count), ],
    3
  )
)

# Create amino-acid composition bar plot
fig1 <- ggplot(
  aa_results,
  aes(AminoAcid, Count)
) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = paste(gene_symbol, "amino-acid composition"),
    x = "Amino acid",
    y = "Number of residues"
  )

print(fig1)

# Save first figure
ggsave(
  "figures/01_amino_acid_composition.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)


# STEP 5. CLINVAR VARIANT DOCUMENTATION
# This step was completed manually using the ClinVar record.
#
# Gene: CSRP3
# Transcript: NM_003476.5
# Variant: c.172T>G (p.Cys58Gly)
# ClinVar Variation ID: 8777
#
# Record the current classification, review status,
# associated conditions, variant type, and publications.
#
# ClinVar link:
# https://www.ncbi.nlm.nih.gov/clinvar/variation/8777/


# STEP 6. SECOND BIOLOGICAL VISUALIZATION
# Plot the mutation position in the reference protein.

if (nzchar(variant_protein)) {
  
  matched <- regmatches(
    variant_protein,
    regexpr("[0-9]+", variant_protein)
  )
  
  position <- as.integer(matched)
  
  if (!is.na(position) &&
      position >= 1 &&
      position <= protein_length) {
    
    fig2 <- ggplot() +
      
      geom_segment(
        aes(
          x = 1,
          xend = protein_length,
          y = 0,
          yend = 0
        ),
        linewidth = 2
      ) +
      
      geom_point(
        aes(x = position, y = 0),
        size = 4
      ) +
      
      annotate(
        "text",
        x = position,
        y = 0.25,
        label = variant_protein,
        size = 4
      ) +
      
      xlim(0, protein_length) +
      ylim(-0.6, 0.6) +
      theme_void() +
      
      labs(
        title = paste(gene_symbol, "reference protein"),
        x = "Amino-acid position"
      ) +
      
      theme(
        axis.title.x = element_text()
      )
    
    print(fig2)
    
    # Save second figure
    ggsave(
      "figures/02_variant_position.png",
      fig2,
      width = 8,
      height = 2.5,
      dpi = 300
    )
  }
}

