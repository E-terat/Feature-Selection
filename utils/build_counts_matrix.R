## dedupmatrix construction
build_counts_matrix <- function(folder,
                                filePattern=".*counts.txt",
                                samplePattern = "PP...")
{
  files <- list.files(path = folder, pattern = filePattern)
  samples <- str_extract(files, pattern = samplePattern)
  
  countsMatrix <- read.delim(paste0(folder,files[1]))
  colnames(countsMatrix) <- c("GeneID", files[1])
  
  for (file in files[2:length(files)]) {
    tmp <- read.delim(paste0(folder,file))
    colnames(tmp) <- c("GeneID", file)
    countsMatrix <- full_join(countsMatrix, tmp, by="GeneID")
  }
  colnames(countsMatrix) <- c("GeneID", samples)
  
  countsMatrix
}
