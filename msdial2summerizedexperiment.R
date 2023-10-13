library(SummarizedExperiment)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)

msdial2se <- function(tablefile){
    tbl <- readr::read_delim(tablefile, delim = "\t", escape_double = FALSE,
                             col_names = FALSE, trim_ws = TRUE)

    numof_lcms_peakmetadata_cols <- 34
    peakmetadata_cols <- tbl %>%
        dplyr::select(1:(numof_peakmetadata_cols-1)) %>%
        dplyr::slice(5:dplyr::n())
    colnames(peakmetadata_cols) <- as.character(unlist(peakmetadata_cols[1,]))
    row_tbl <- peakmetadata_cols[-1,]

    quantval_cols <- tbl %>% dplyr::select(numof_peakmetadata_cols:ncol(tbl))

    numof_sample_class <- length(na.omit(unique(unlist(quantval_cols[1,]))))
    numof_avgstdev_cols <- numof_sample_class * 2
    quantval_cols <- quantval_cols %>%
        dplyr::select(1:(ncol(quantval_cols)-numof_avgstdev_cols))

    sample_meta_tbl <- t(quantval_cols[c(1,2,3,4,5),])
    colnames(sample_meta_tbl) <- c("Class", "FileType",
                                   "InjectionOrder", "BatchId", "SampleId")
    col_tbl <- tibble::as_tibble(sample_meta_tbl)
    #col_data <- tibble::column_to_rownames(col_tbl, var = "SampleId")

    quantval_tbl <- quantval_cols[-c(1,2,3,4,5),] %>%
        dplyr::mutate(across(where(is.character), as.numeric))

    se <- SummarizedExperiment::SummarizedExperiment(assays = quantval_tbl,
                                                     rowData = row_tbl,
                                                     colData = col_tbl)
    return(se)
}
