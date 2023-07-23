# functions used in the analysis
filter_expression <- function(countMatrix,
                             CPMcutoff=0.5,
                             minSamples=0.75)
{
  keep <- rowSums(edgeR::cpm(countMatrix) > CPMcutoff) >= (floor(ncol(countMatrix)*minSamples))
  keep
}

norm_vst <- function(countData,
                     colData,
                     formula= ~ 1){
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = colData,
                                        formula)
  vsd <- DESeq2::vst(dds, blind = FALSE)
  vsd
}