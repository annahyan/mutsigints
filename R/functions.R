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
    
    odds.mat [ p.out > p.val.threshold ] = 0
    
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
#' @param robust Specifies if robust statistics should be used. 
#' Robust functions are called from robustbase package. Default: TRUE
#' @param path.to.sig If true, linear regression is used to predict signatures
#' from pathways. Otherwise, logistic regression is used to predict pathways 
#' from signatures. Default: TRUE
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @return A matrix with signature-number of rows and pathway-number of columns
#' with values indicating odds ratio of corresponding signatures and pathways.
#' 

get_sig_path_lms = function(sigs.df, pathways.df, 
                            sig.log = TRUE, 
                            robust = TRUE,
                            path.to.sig = TRUE,
                            p.val.threshold = 0.05, 
                            p.adjust = TRUE, method = "BH",
                            ...) {
    
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
    p.values = matrix(1, nrow = length(sigs), ncol = length(pathways),
                      dimnames = list(sigs, pathways))
    
    for (sig in sigs) {
        for (pathway in pathways) {
            # cat(sig, pathway, "\n")
            
            ## Applying a threshold on minimal number of non-zero elements
            # cat("before printing zero element counts.\n")
            zero.sigs = sum(tissue.concat[, sig] != 0)
            zero.paths = sum(tissue.concat[, pathway] != 0)
            
            if (zero.sigs < 3 | zero.paths < 3) {
                # cat("zero.sigs = ", zero.sigs, " zero paths = ", zero.paths, "\n")
                int.mat[sig, pathway] = 0
                p.values[sig, pathway] = 1
                next
            }
            
            paths.binary = as.numeric(tissue.concat[, pathway] > 0)
            
            ### Drop if all elements are non-zero
            if (length(unique(paths.binary)) == 1) {
                int.mat[sig, pathway] = 0
                p.values[sig, pathway] = 1
                next
            }
            
            if ( robust ) {
                if (path.to.sig) {
                    ## with Robustbase
                    rob.lin.mod = robustbase::lmrob(
                        tissue.concat[, sig] ~ tissue.concat[, pathway], ...)
                    int.mat[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Estimate"][2]
                    p.values[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Pr(>|t|)"][2]
                    ## with rlm
                    # rob.lin.mod = MASS::rlm(tissue.concat[, sig] ~ tissue.concat[, pathway], ...)
                    # int.mat[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Value"][2]
                    # p.values[sig, pathway] = tryCatch({
                    #     sfsmisc::f.robftest(rob.lin.mod, var = -1)$p.value},
                    #     error = function(e) {return(1)})
                    
                } else {
                    ### Trying robust logistic regression first.
                    ### If it throws an error, then regular logistic regression.
                    
                    log.mod = tryCatch({
                        rob.out = robustbase::glmrob(
                            paths.binary ~ 1 + tissue.concat[, sig], 
                            family = binomial)
                        rob.out
                        }, 
                        error = function(e) {
                            print(e$message)
                            cat("Robust regression cannot be computed. \n Running regular regression instead.\n")
                            glm.out = glm(
                                paths.binary ~ 1 + tissue.concat[, sig], 
                                family = binomial)
                            glm.out
                        })
                    int.mat[sig, pathway] = summary(log.mod)$coefficients[, "Estimate"][2]
                    p.values[sig, pathway] = summary(log.mod)$coefficients[, "Pr(>|z|)"][2]
                }
            } else {
                if (path.to.sig) {
                    # lin.mod = lm(tissue.concat[, sig] ~ tissue.concat[, pathway])
                    lin.mod = lm(tissue.concat[, sig] ~ paths.binary)
                    int.mat[sig, pathway] = summary(lin.mod)$coefficients[, "Estimate"][2]
                    p.values[sig, pathway] = summary(lin.mod)$coefficients[,"Pr(>|t|)"][2]
                } else {
                    paths.binary = as.numeric(tissue.concat[, pathway] > 0)
                    log.mod = glm(
                        paths.binary ~ 1 + tissue.concat[, sig], 
                        family = binomial, ...)
                    int.mat[sig, pathway] = summary(log.mod)$coefficients[, "Estimate"][2]
                    p.values[sig, pathway] = summary(log.mod)$coefficients[, "Pr(>|z|)"][2]
                }
            }
        }
    }
    
    p.out = p_matrix_process(p.values, p.adjust = p.adjust, method = method)

    int.mat[p.out > p.val.threshold] = 0
    # return(list(int.mat = int.mat, p.value = p.out))
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



#' Assess signature-pathway interactions across tissues for PCAWG with a custom
#' function
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
        
        tissue.odds.plot = myggheatmap(
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



#' For a list of interaction metrics this function summarizes the positive and
#' negative interactions in a plot.
#' 
#' @param list.of.int.elems A list with matrix elements.
#' @param threshold All the values below the threshold are discarded.
#' 
#' @details The input is a list of matrices corresponding to e.g. interactions
#' in different tissues. The output plot summarizes the number of list elements
#' in which this interaction happens either in positive or negative direction.
#' The positive and negative interactions are summarized separately.
#' 
#' @return ggplot object

plot_all_counts = function(list.of.int.elems, threshold = 0.1) {
    
    # all.sigs = do.call(c, lapply(summary.matrix, function(x) colnames(x)) ) %>%
    #     unique()
    
    all.rows = do.call(c, lapply(list.of.int.elems, function(x) rownames(x))) %>% 
        unique()
    
    all.cols = do.call(c, lapply(list.of.int.elems, function(x) colnames(x))) %>% 
        unique()
    
    pos.ints.mat = matrix(0, nrow = length(all.rows), ncol = length(all.cols),
                          dimnames = list(all.rows, all.cols))
    
    neg.ints.mat = matrix(0, nrow = length(all.rows), ncol = length(all.cols),
                          dimnames = list(all.rows, all.cols))
    
    for (rowelem in all.rows) {
        for (colelem in all.cols) {
            point.ints = unlist(sapply(list.of.int.elems, 
                                       function(x) as.data.frame(x)[rowelem, colelem]))
            
            point.ints[abs(point.ints) < threshold] = 0
            
            pos.ints = sum(point.ints > 0, na.rm = TRUE)
            neg.ints = -sum(point.ints < 0, na.rm = TRUE)
            
            pos.ints.mat[rowelem, colelem] = pos.ints
            neg.ints.mat[rowelem, colelem] = neg.ints
        }
    }
    
    pos.pivot = pos.ints.mat %>% 
        as.data.frame() %>% rownames_to_column(var = "rows") %>% 
        pivot_longer(cols = -rows, values_to = "count", names_to = "cols") %>% 
        mutate(int.type = "pos")
    
    # rows   name       count int.type
    # <chr>  <chr>      <dbl> <chr>   
    # 1 Ageing Cell Cycle     1 pos     
    # 2 Ageing HIPPO          2 pos     
    # 3 Ageing NRF2           1 pos     
    # 4 Ageing PI3K           2 pos   
    
    neg.pivot = neg.ints.mat %>% 
        as.data.frame() %>% rownames_to_column(var = "rows") %>% 
        pivot_longer(cols = -rows, values_to = "count", names_to = "cols") %>% 
        mutate(int.type = "neg")
    
    
    gg.final.dt = bind_rows(pos.pivot, neg.pivot)
    
    gg.final.dt = as.data.frame(gg.final.dt)
    
    abs.row.nonzero = gg.final.dt %>% 
        group_by(rows) %>% 
        summarise(abssum = sum(abs(count))) %>% 
        filter(abssum > 0) %>% 
        pull(rows)
    
    abs.col.nonzero = gg.final.dt %>% 
        group_by(cols) %>% 
        summarise(abssum = sum(abs(count))) %>% 
        filter(abssum > 0) %>% 
        pull(cols)
    
    
    row.indices = setNames(1:length(abs.row.nonzero), abs.row.nonzero)
    col.indices = setNames(1:length(abs.col.nonzero), abs.col.nonzero)
    
    gg.final.dt = gg.final.dt %>% 
        filter(rows %in% abs.row.nonzero) %>% 
        filter(cols %in% abs.col.nonzero) %>% 
        mutate(x = row.indices[rows],
               y = col.indices[cols])
    
    
    gg.final.dt = gg.final.dt %>%
        mutate(xlab = ifelse(int.type == "pos", x - 0.2, x + 0.2),
               ylab = ifelse(int.type == "pos", y + 0.15, y - 0.15),
               int.type = factor(int.type, levels = c("pos", "neg") ) )

    # gg.final.dt = gg.final.dt %>% 
    #     mutate(xlab = x, 
    #            ylab = y)
    
    gg.final.dt = as.data.frame(gg.final.dt)
    
    d <- ggplot(gg.final.dt, aes(x = x, y = y, color = int.type,
             shape = int.type, label = count, alpha = abs(count)))
    
    if (Sys.info()['sysname'] == "Darwin") {
        # d = d + point_with_family(geom_point(size = 5), "wqy-microhei")
        d = d + point_with_family(geom_point(size = 1), "Arial Unicode MS")
    } else {
        d = d + geom_point(size = 4.5)
    }
    d = d +
        point_with_family(geom_point(size = 18), "impact") +
        geom_text(aes(x = xlab, y = ylab), size = 2, color = "black", fontface = "bold") +
        scale_shape_manual(values=c("\u25E4","\u25E2")) +
        scale_color_brewer(palette = "Set1") 
    
    d = d +
        # theme_minimal() +
        theme( panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.ticks = element_blank(),
               # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
               axis.text.x = element_text(
                   angle = 45, 
                   margin = margin(r = 0, t = 0, b = 0, l = 0),
                   hjust = 0),
               axis.text.y = element_text(
                   margin = margin(r = 0, t = 0, b = 0, l = 0)),
               axis.title = element_blank(),
               legend.position = "none"
        ) +
        scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)),
            breaks = row.indices,
            labels = names(row.indices),
            position = "top") +
        scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),
            breaks = col.indices,
            labels = names(col.indices))
    d
    return(d)
}



#' Summarizing individual matrices of interactions
#' 
#' @param interaction.list Input list
#' @param pos.ints Defines if positive or negative interactions 
#' should be considered. Default: TRUE
#' @return Matrix with element is concatenated tissue names
#' bearing the interaction 
#' 
summarize_int_mat = function(interaction.list, pos.ints = TRUE) {
    
    interaction.list = lapply(interaction.list, as.data.frame)
    
    active.colnames = do.call(c, sapply(interaction.list, colnames) ) %>%
        unique() %>% sort
    active.rownames = do.call(c, sapply(interaction.list, rownames) ) %>%
        unique() %>% sort
    
    outmat = matrix("", nrow = length(active.rownames), 
                    ncol = length(active.colnames),
                    dimnames = list(active.rownames,
                                    active.colnames))
    
    for(rowelem in active.rownames) {
        for(colelem in active.colnames) {
            vals = sapply(interaction.list, function(x) 
                as.data.frame(x)[rowelem, colelem]) %>% unlist
            if(pos.ints) {
                mat.names = which(vals > 0) %>% names
            } else {
                mat.names = which(vals < 0) %>% names
            }
            outmat [rowelem, colelem] = paste(mat.names, collapse = ",")
        }
    }
    return(outmat)
}

### https://stackoverflow.com/questions/48531257/ggplot-geom-point-how-to-set-font-of-custom-plotting-symbols
point_with_family <- function(layer, family) {
    old_geom <- layer$geom
    new_geom <- ggproto(
        NULL, old_geom,
        draw_panel = function(self, data, panel_params, coord, na.rm = FALSE) {
            pts <- ggproto_parent(GeomPoint, self)$draw_panel(
                data, panel_params, coord, na.rm = na.rm
            )
            pts$gp$fontfamily <- family
            pts
        },
        draw_key = function(self, data, params, size) {
            pts <- ggproto_parent(GeomPoint, self)$draw_key(
                data, params, size
            )
            pts$gp$fontfamily <- family
            pts
        }
    )
    layer$geom <- new_geom
    layer
}
