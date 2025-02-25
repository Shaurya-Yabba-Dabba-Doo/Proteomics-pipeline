library(RDAVIDWebService)

david <- DAVIDWebService$new(email = "your_email@example.com")

protein_list <- read.table("your_protein_list.txt", header = TRUE, stringsAsFactors = FALSE)

protein_vector <- as.character(protein_list$ProteinID)

addList(david, protein_vector, idType = "UNIPROT_ACCESSION", listName = "MyProteinList")

annotations <- getFunctionalAnnotation(david)

enrichment_results <- getFunctionalAnnotationChart(david)

write.table(annotations, "david_annotations.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(enrichment_results, "david_enrichment_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

removeList(david, "MyProteinList")

close(david)
