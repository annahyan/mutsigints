#' The function aggregates mutations based on annotations
#' @param annotations annotations data frame. It contains must contain 2 columns-
#' Signature and Annotation
#' @param signature_matrix_df Signature matrix, which should be annotated
#' @return A matrix with same number of rows as the input siganture_matrix_df and 
#' aggregated and renamed columns according to annotations dataframe.


set_signature_annotation = function(annotations, signature_matrix_df) {
    
    sig_annotation_groups = annotations$Annotation [
        match(colnames(signature_matrix_df), 
              annotations$Signature) ]
    
    for(i in 1:length(sig_annotation_groups)) {
        if (is.na(sig_annotation_groups[i]) ) {
            sig_annotation_groups[i] = colnames(signature_matrix_df)[i]
        }
    }
    
    signature_matrix_df.ann = signature_matrix_df %>% t() %>% 
        rowsum(group = sig_annotation_groups) %>% t()
    
    return(signature_matrix_df.ann)
}
