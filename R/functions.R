#' The function aggregates mutations based on annotations
#' @param annotations annotations data frame. It contains must contain 2 columns-
#' Signature and Annotation
#' @param signature.matrix.df Signature matrix, which should be annotated
#' @return A matrix with the same number of rows as the input 
#' siganture.matrix.df and aggregated and renamed columns 
#' according to annotations dataframe.


set_signature_annotation = function(annotations, signature.matrix.df) {
    
    sig.annotation.groups = annotations$Annotation [
        match(colnames(signature.matrix.df), 
              annotations$Signature) ]
    
    for(i in 1:length(sig.annotation.groups)) {
        if (is.na(sig.annotation.groups[i]) ) {
            sig.annotation.groups[i] = colnames(signature.matrix.df)[i]
        }
    }
    
    signature.matrix.df.ann = signature.matrix.df %>% t() %>% 
        rowsum(group = sig.annotation.groups) %>% t()
    
    return(signature.matrix.df.ann)
}


#' Performs p-value adjustment on a matrix of p-values based on input parameters. 
#' 
#' @return A matrix with the same dimensions and dimnames as the input
#' 
p_matrix_process = function(p.values, p.adjust = TRUE, method = "BH") {
    if (p.adjust) {
        p.out <- p.values %>% 
            as.vector %>% p.adjust(method = method) %>% 
            matrix(ncol = ncol(p.values), dimnames = dimnames(p.values))
    } else {
        p.out = p.values
    }
    
}




#' Remove the columns and rows of the dataframe with all 0's
#' 
rm_zeros = function(x) {
    pos.x = abs(x)
    rmout = x[rowSums(pos.x) > 0, colSums(pos.x) > 0, drop = FALSE]
    return(rmout)
}



#' Extract tissue and the signatures active in that tissue from a matrix of all 
#' the data.
#' 
#' @param sigs.full.subset is the signature activity matrix, where first column
#' is the cancer type/tissue and signature activities start from column 4. This
#' is the input type for signature activities from Alexandrov et al.(Nature, 2020)
#' paper.
#' @param tissue Tissue of interest
#' 
#' @return Returns a subset of samples from the selected tissue and only active 
#' signatures
#' 
subset_tissue = function(sigs.full.subset, tissue) {
    
    tissue.signatures = sigs.full.subset[ sigs.full.subset[, 1] == tissue, ]
    
    signature.matrix = tissue.signatures[, 4:ncol(sigs.full.subset)]
    signature.matrix = signature.matrix[, colSums(signature.matrix) > 0, 
                                        drop = FALSE]
    out = cbind(tissue.signatures[, 1:3], signature.matrix) 
    return(out)
}


#' For a given signature activity matrix and pathway alteration matrix 
#' computes Fisher's exact test for all pathway-signature interactions.
#' @param sigs.df signatures dataframe or matrix
#' @param pathways.df pathways dataframe or matrix
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @return A matrix with signature-number of rows and pathway-number of columns
#' with values indicating odds ratio of corresponding signatures and pathways.
#' 
get_sig_path_assocs = function(sigs.df, pathways.df, p.val.threshold = 0.05, 
                               p.adjust = TRUE, method = "BH") {
    
    sigs = colnames(sigs.df)
    pathways = colnames(pathways.df)
    
    if (! all(rownames(sigs) == rownames(pathways))) {
        stop("Colnames corresponding to sample names should match between two matrices.")
    }
    
    odds.mat = matrix(0, nrow = length(sigs), ncol = length(pathways),
                      dimnames = list(sigs, pathways))
    
    p.values = matrix(0, nrow = length(sigs), ncol = length(pathways),
                      dimnames = list(sigs, pathways))
    
    for (sig in sigs) {
        for (pathway in pathways) {
            
            dd = cbind(sigs.df[, sig, drop = FALSE],
                       pathways.df[, pathway, drop = FALSE])
            dd = na.omit(dd)
            
            dd[dd > 0] = 1
            
            dd[, pathway] = factor(dd[ ,pathway], levels = c(0, 1), 
                                   labels = c("WT", "Mutated"))
            dd[ , sig] = dd[, sig] > 0
            dd[, sig] = factor(dd[, sig], levels = c(FALSE,TRUE), 
                               labels = c("Inactive", "Active"))
            fisher.out = fisher.test(table(dd) + 1) ### adding 1, to avoind +/- Inf values
            odds.mat[sig, pathway] = log2(fisher.out$estimate)
            p.values[sig, pathway] = fisher.out$p.value
        }
    }
    
    p.out = p_matrix_process(p.values, p.adjust = p.adjust, method = method)
    
    # pheatmap(p.out)
    
    odds.mat [ p.out >= p.val.threshold ] = 0
    
    # odds.mat = odds.mat[rowSums(abs(odds.mat), na.rm = TRUE) > 0, , drop = FALSE]
    # odds.mat = odds.mat[, colSums(abs(odds.mat), na.rm = TRUE) > 0 , drop = FALSE]
    
    return(as.data.frame(odds.mat))
}



#' For signature and pathway matrices performs linear regression for 
#' each signature-pathway pair.
#' For a given signature activity matrix and pathway alteration matrix 
#' computes Fisher's exact test for all pathway-signature interactions.
#' @param sigs.df signatures dataframe or matrix
#' @param pathways.df pathways dataframe or matrix
#' @param sig.log Defines whether signatures should be logged. Default: TRUE
#' @param robust Specifies if robust statistics should be used. Default: TRUE
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @return A matrix with signature-number of rows and pathway-number of columns
#' with values indicating odds ratio of corresponding signatures and pathways.
#' 

get_sig_path_lms = function(sigs.df, pathways.df, 
                            sig.log = TRUE, 
                            robust = TRUE,
                            p.val.threshold = 0.05, 
                            p.adjust = TRUE, method = "BH") {
    
    if (sig.log) {
        sigs.df = sigs.df %>% 
            mutate(across(.cols = everything(), ~ log(.x + 1 )))
    }
    
    sigs = colnames(sigs.df)
    pathways = colnames(pathways.df)
    
    if( ! all.equal(rownames(sigs), rownames(colnames)) ) {
        stop("The rownames (corresponding to sample names) in both datasets should match.")
    }
    
    tissue.concat = merge(pathways.df, sigs.df, 
                          by = "row.names")
    
    int.mat = matrix(0, nrow = length(sigs), ncol = length(pathways),
                     dimnames = list(sigs, pathways))
    p.values = matrix(0, nrow = length(sigs), ncol = length(pathways),
                      dimnames = list(sigs, pathways))
    
    for (sig in sigs) {
        for (pathway in pathways) {
            # cat(sig, pathway, "\n")
            
            if ( robust ) {
                rob.lin.mod = MASS::rlm(tissue.concat[, sig] ~ tissue.concat[, pathway])
                int.mat[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Value"][2]
                p.values[sig, pathway] = tryCatch({
                    sfsmisc::f.robftest(rob.lin.mod, var = -1)$p.value},
                    error = function(e) {return(1)})
            } else {
                lin.mod = lm(tissue.concat[, sig] ~ tissue.concat[, pathway])
                int.mat[sig, pathway] = summary(lin.mod)$coefficients[, "Estimate"][2]
                p.values[sig, pathway] = summary(lin.mod)$coefficients[,"Pr(>|t|)"][2]
            }
        }
    }
    
    p.out = p_matrix_process(p.values, p.adjust = p.adjust, method = method)

    int.mat[p.out < p.val.threshold] = 0
    return(int.mat)
}

#' Extract signature and pathway activities for a tissue
#' @param tissue Tissue to extract
#' @param sigs.input Signature activities
#' @param pathways.input Pathway activities 
#' @return A list with sigs and pathways elements with corresponding matrices.
#' Can be forwarded do interaction calculting functions

get_tissue_pathway_activities = function(tissue, 
                                         sigs.input, 
                                         pathways.input) {
    
    tissue.sig.subset = subset_tissue(sigs.input, tissue = tissue) %>% 
        select(4:ncol(.))
    
    tissue.rownames = rownames(tissue.sig.subset)  
    
    tissue.path.subset = pathways.input[ tissue.rownames, ] 
    tissue.path.subset = tissue.path.subset %>% 
        select(4:ncol(.)) %>% 
        select_if(colSums(., na.rm = TRUE) != 0)
    
    return(list(sigs = tissue.sig.subset, paths = tissue.path.subset))
}



#' Assess signature-pathway interactions across tissues for PCAWG
#' @param sigs.input Signature activities
#' @param pathways.input Pathway status - formatted like mutated.pathways.tissues
#' Pathway activities start at column 4.
#' @param interaction_function The function defining the metric. E.g. get_sig_path_assocs
#' @param path.min.tissues Minimal number of samples in each tissue to be considered 
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @param plot.limits Scale limits in the ggheatmap
#' @return A list with the length of abundant tissues in the datasets, 
#' where each element is the interaction matrix

pcawg_sig_pathway_int = function(sigs.input,
                                 pathways.input,
                                 interaction_function, 
                                 path.min.tissues = 30,
                                 p.val.threshold = 0.1,
                                 p.adjust = TRUE,
                                 method = "BH",
                                 ...) {
    
    abundant.tissues = which(sigs.input$Cancer.Types %>% 
                                 table() > path.min.tissues) %>% names()
    
    tissue.odds.mats = list()
    
    for (tissue in abundant.tissues){
        cat(tissue, "\n")
        
        tissue.elems = get_tissue_pathway_activities(tissue, 
                                                     sigs.input, 
                                                     pathways.input)

        tissue.odds.mat = interaction_function(
            tissue.elems$sigs, 
            tissue.elems$paths,
            p.val.threshold = p.val.threshold, 
            p.adjust = p.adjust,
            method = method, ...) %>% 
            rm_zeros()
        
        tissue.odds.mats[[tissue]] = tissue.odds.mat
    }
    return(tissue.odds.mats)
}

# A ggheatmap wrapper for interaction matrix returns
#' @param metric.matrix the input matrix
#' @description All other parameters are passed to ggheatmap
#' \link[ggheatmap]{ggheatmap} 
#' @return Returns a ggheatmap (which is a gg object, 
#' but not a conventional gg object unfortunately.)

ggheatmap_wrapper = function(metric.matrix,
                             colorPalette = "RdBu",
                             points = T, revColors = F, revRow = T, 
                             scaleName = "log(OR)",
                             title = NULL,
                             limits = NULL, 
                             ...) {
    
    if (( ! all(dim(metric.matrix) > 0) ) | is.null(dim(metric.matrix) ) ) {
        tissue.odds.plot = NA
        cat("The input matrix has no active dimensions or is null.")
    } else {
        orderrow = T
        ordercol = T
        if(dim(metric.matrix)[1] == 1) {
            orderrow = F
        } 
        
        if(dim(metric.matrix)[2] == 1){
            ordercol = F
        } 
        
        if (is.null(limits)) {
            range_lim = round(max(abs(range(metric.matrix))) * 1.1, 1)
            limits = c(-range_lim, range_lim)
        }
        
        tissue.odds.plot = ggheatmap::ggheatmap(
            as.data.frame(metric.matrix), colorPalette = colorPalette,
            orderCol = ordercol,
            orderRow = orderrow,
            points = points, 
            revColors = revColors, 
            revRow = revRow, 
            scaleName = scaleName,
            title = title,
            limits = limits, ...)
    }
    
    return(tissue.odds.plot)
}
