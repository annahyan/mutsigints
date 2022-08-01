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
#' @param psize Controls the size of the triangles. Default: 8.
#' @param lsize Label size. Default: 2.
#' @param expand.mult A vector of expansion factors for x and y axes. 
#' Default. c(0.04, 0.04)
#' @param diag.only If true, for symmetric matrices only diagonal will be printed.
#' Default:TRUE.
#' 
#' @details The input is a list of matrices corresponding to e.g. interactions
#' in different tissues. The output plot summarizes the number of list elements
#' in which this interaction happens either in positive or negative direction.
#' The positive and negative interactions are summarized separately.
#' 
#' @return ggplot object

plot_all_counts = function(list.of.int.elems, threshold = 0.1, 
                           psize = 8, lsize = 2, expand.mult = c(0.04, 0.04), 
                           diag.only = TRUE) {
    
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
    
    ### If the matrices are symmetric, then only upper triangle is plotted.

    if (all.equal(neg.ints.mat, t(neg.ints.mat)) &
            all.equal(pos.ints.mat, t(pos.ints.mat)) & diag.only) {
        gg.final.dt = gg.final.dt %>% filter(x < y)
        
        row.indices = row.indices[ which(row.indices %in% gg.final.dt$x)]
        col.indices = col.indices[ which(col.indices %in% gg.final.dt$y)]
    }
    
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
        d = d + point_with_family(geom_point(size = psize), "Arial Unicode MS")
    } else {
        d = d + geom_point(size = psize)
    }
    d = d +
        # point_with_family(geom_point(size = 18), "impact") +
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
        scale_x_continuous(expand = expansion(mult = expand.mult),
            breaks = row.indices,
            labels = names(row.indices),
            position = "top") +
        scale_y_continuous(expand = expansion(mult = expand.mult),
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


#' Calculates interaction values using a given metric_func in a dataset
#' through bootstrapping.
#' @param data.input Input mutational signatures
#' @param metric_func Function calculating the interaction from a matrix
#' @param min.tissue.samples Tissues with less than 
#' this number of samples are discarded. Default: 20
#' @param sample.rate The proportion of samples which should be sampled during 
#' the bootstrapping procedure. Default: 0.9
#' @param sample.N Number of times to sample. Default: 100
#' @param N Number of times the function should be called for each sampled set:
#' Default: 1
#' @param seed Random number generator seed. Default: 1.
#' @param ... Parameters to be passed to metric_func.

get_metrics_list = function(data.input, metric_func,
                            min.tissue.samples = 20,
                            sample.rate = 0.9, sample.N = 100, 
                            N = 1, seed = 1,  ...) {
    
    action_func = function(x) {
        metric_func(x, ...)
    }
    
    tissues = unique(data.input$Cancer.Types)
    
    int.metric.list = list()
    
    for (tissue in tissues) {

        #### Filtering for tissues with enough samples
        tissue.signatures = list()
        cat("Processing tissue ", tissue, "\n")
            
        dt = data.input
        tissue.sig = dt %>% filter(Cancer.Types == tissue) %>%
            select(4:ncol(dt))
        # tissue.sig = tissue.sig %>% filter(MMR == 0)
        tissue.sig = tissue.sig[, colSums(tissue.sig) > 10,
                                drop = FALSE]
        
            if (nrow(tissue.sig) < min.tissue.samples) {
                cat("Tissue:", tissue , " has <", min.tissue.samples,
                    "samples in ", dataset, ". Skipping.\n")
                tissue.signatures = NULL
            } else {
                tissue.signatures = tissue.sig
            }
        
        metric_out = run_calc_metrics(tissue.signatures, action_func,
                                      sample.rate, sample.N, N, seed)
        int.metric.list[[tissue]] = metric_out
    }
    return(int.metric.list)
}


# Calculates the metric for a list or a dataframe.

run_calc_metrics = function(tissue.signatures, action_func, sample.rate = 0.9, 
                            sample.N = 100, N = 1, seed = 1) {
    
    if (inherits(tissue.signatures, "list")) {
        int.metric.calcs = lapply(tissue.signatures, function(dt)
            calc_metric(dt, action_func, sample.rate, sample.N, N, seed)
        )
    } else if (inherits(tissue.signatures, "data.frame")) {
        int.metric.calcs = calc_metric(tissue.signatures, action_func, 
                                       sample.rate, sample.N, N, seed)
    }
    return(int.metric.calcs)
}

#' Calculates the metric function for a matrix, by performing bootstrapping.
#' @param dt.sigs A dataframe or a matrix of mutational signatures.
#' @param action_func Function calculating the interaction from a matrix
#' @param sample.rate The proportion of samples which should be sampled during 
#' the bootstrapping procedure. Default: 0.9
#' @param sample.N Number of times to sample. Default: 100
#' @param N Number of times the function should be called for each sampled set:
#' Default: 1
#' @param seed Random number generator seed. Default: NULL.

calc_metric = function(dt.sigs, action_func, sample.rate, sample.N, N, seed = NULL ) {
    
    if (!is.null(seed)) { set.seed(seed)}
    if (is.null(dt.sigs)) {
        return(NULL)
    } else {
        sample.length = nrow(dt.sigs)
        sample.counts = round(sample.length * sample.rate)
        
        sampled.out = lapply(1:sample.N, function(i) {
            sampled.signatures = sample_n(dt.sigs, sample.counts)
            ll <- replicate(N, action_func(sampled.signatures),
                            simplify = FALSE )
            
            ll.named = lapply(ll, function(x) {
                colnames(x) = colnames(sampled.signatures)
                rownames(x) = colnames(sampled.signatures)
                return(x)
            })
            return(ll.named)
            
        }) }
}

#' Extracts the interaction metrics of a tissue 
#' from a list of interaction networks.
#' @param tissue Tissue to be extracted
#' @param network.lists The input list
#' @param filter.list A list of metric values to be filtered out. The list 
#' element names have the same names as metric names in network.lists. 
#' The values specify filtering threshold. Everything below the number is 
#' set to 0. Default: NULL
#' @return A list with interaction metric matrices. The list element names match
#' those in the input network.lists.

get_tissue_dataset_networks = function(tissue,
                                       network.lists, 
                                       filter.list = NULL) {
    
    out.list = lapply(network.lists, function(x) 
        x[[tissue]])
    
    if (is.null(filter.list)) {
        return(out.list)
    } 
    
    filter.list.absent = setdiff(names(filter.list), names(out.list))
    
    if (length( filter.list.absent > 0))  {
        stop(paste("filter.llist has variable names not present:", 
                   paste(filter.list.absent, collapse = " ")) )
    } 
    
    for (var in names(filter.list)) {
        varlim = filter.list[[var]]
        mat = out.list[[var]]
        mat[ abs(mat) < varlim ] = 0
        mat = mat[rowSums(mat) > 0, colSums(mat) > 0]
        out.list[[var]] = mat
    }
    return(out.list)
}


#' Plotting a heatmap for signature networks. Can be a PCAWG-format with first 
#' 3 columns being a metadata, or a matrix/data.frame with all signature activity
#' values.
#' @param .dataset Input dataset
#' @param filename If provided the pheatmap will be save with this filename.
#' Default: NULL
#' @param main Title of the heatmap
#' @param ... Parameters will be passed to pheatmap function.
#' 
#' 

plot_tissue_heatmap = function(.dataset, filename = NULL, main = NULL, ...) {
    
    classes = sapply(.dataset, class)
    if (! all(sapply(.dataset, class)[1:3] %in% c("numeric", "integer")))  {
        dt.plot = .dataset[, 4:ncol(.dataset)]
    } else {
        
        dt.plot = .dataset
    }
    
    dt.plot = dt.plot[, colSums(dt.plot) > 0]
    
    if (is.null(main)) {
        main = paste(nrow(.dataset), "samples")
    } else {
        main = paste(main, "-", nrow(.dataset), "samples")
    }
    
    if (is.null(filename) ) {
        pheatmap(log(dt.plot + 1), 
                 color = viridis(15),
                 show_rownames = FALSE,
                 main = main, 
                 width = 7.1, height = 6.83,
                 # filename = here("figures/signature_heatmaps/", "TCGA.annotated.tissues.png"),
                 # annotation_row = tcga.annotation.row,
                 # annotation_colors = tcga.mycolors
                 # # annotation_legend = FALSE
                 ...)
    } else {
        pheatmap(log(dt.plot + 1), 
                 color = viridis(15),
                 show_rownames = FALSE,
                 width = 7.1, height = 6.83,
                 main = main,
                 filename = filename,
                 # annotation_row = tcga.annotation.row,
                 # annotation_colors = tcga.mycolors
                 # # annotation_legend = FALSE
                 ...)   
    }
}


#' Adds a legend title to pheatmap legend
#' 
#' @param p ggplot object produced by pheatmap.
#' @param legend.grob.index The grob which contains the legend grob.
#' @details Based on Mike H.'s answer from https://stackoverflow.com/questions/36852101/r-legend-title-or-units-when-using-pheatmap
#' 

add_pheatmap_legend_title = function(p, legend.grob.index = 6, 
                                     legend.text = "log(n)") {
    
    legend.grob <- p$gtable$grob[[legend.grob.index]]
    
    legend.grob$children[[1]]$y <- legend.grob$children[[1]]$y - 0.08 * legend.grob$children[[1]]$y 
    legend.grob$children[[2]]$y <- legend.grob$children[[2]]$y - 0.08 * legend.grob$children[[2]]$y
    legend.grob$children[[1]]$x <- legend.grob$children[[1]]$x + 0.15 * legend.grob$children[[1]]$x
    legend.grob$children[[2]]$x <- legend.grob$children[[2]]$x + 0.15 * legend.grob$children[[2]]$x
    
    leg_label <- textGrob(legend.text,
                          x=0,y=0.95,hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))
    
    legend.grob2 <- addGrob(legend.grob,leg_label)
    
    p$gtable$grobs[[legend.grob.index]] = legend.grob2
    return(p)
}


#' The function summarizes interactions for a given metric across tissues, 
#' identifies common interaction motifs.
#' @param metric.list input metric.list, element names are different interaction
#' metrics. Each element is a list corresponding to interactions corresponding
#' to individual tissues.
#' @param metric Metric for which the summaries should be obtained.
#' @param outdir The output directory where the interaction matrices and the 
#' outputs should be written.
#' @param threshold All interactions below this threshold are set to 0. 
#' Default: 0.2.
#' 

get_common_sigs = function(metric.list, metric, outdir, threshold = 0.2) {

    for (tissue in names(metric.list[[metric]])) {
        cat("\t\tTissue = ", tissue, "\n")
        tissue.int.mat = metric.list[[metric]][[tissue]]  
        
        if(is.null(tissue.int.mat)) next
        
        tissue.int.mat[ abs(tissue.int.mat) < threshold] = 0
        tissue.int.mat  = sign(tissue.int.mat)
        write.table(tissue.int.mat, file = paste0(outdir, "/",  tissue, ".txt"),
                    sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
    }
    python.call = paste("python3", here("python", "combine_networks.py"),
                        outdir)
    cat(python.call, "\n")
    system(python.call, wait = TRUE)
}


#' 
#' @param network.lists The input list.
#' @param tissue Tissue to be extracted.
#' @param filter.list A list of metric values to be filtered out. The list 
#' element names have the same names as metric names in network.lists. 
#' The values specify filtering threshold. Everything below the number is 
#' set to 0. Default: NULL

concat_networks = function(network.lists, tissue, filter.list) {
    
    try({
        tissue.nets = get_tissue_dataset_networks(
            tissue,
            network.lists = network.lists, 
            filter.list = filter.list
        )
        
        for (type in names(tissue.nets)) {
            if (sum(abs(tissue.nets[[type]]) ) == 0  ) {
                tissue.nets[[type]] = NULL
            }
        }
        
        ts.multi.graph = tissue_multi_graph(tissue.nets)
        ## layout = "dh"
        ts.adj.mat = as.matrix(as_adjacency_matrix(ts.multi.graph))
        ts.adj.mat[ts.adj.mat < 2] = 0
        ts.adj.mat[ts.adj.mat > 0] = 1
        ts.adj.mat = ts.adj.mat + t(ts.adj.mat)
        
        
        edge.info = ts.multi.graph %>% activate(edges) %>% data.frame() %>% 
            filter(int_type != "MI")
        
        adj.net = as.matrix(as_adjacency_matrix(ts.multi.graph)) * 0
        for(i in 1:nrow(edge.info)) {
            adj.net[edge.info[i, "from"], edge.info[i, "to"]] = 
                adj.net[edge.info[i, "from"], edge.info[i, "to"]] + edge.info[i, "weight_"]
        }
        adj.net = adj.net + t(adj.net)
        ts.adj.mat = ts.adj.mat * adj.net
        ts.adj.mat = ts.adj.mat [colSums(abs(ts.adj.mat) ) > 0, 
                                 rowSums(abs(ts.adj.mat)) > 0]
        return(ts.adj.mat)
    })
    
    return(c())
}



#' Wrapper of tissue interaction network plot
#' @param tissue Tissue to be extracted
#' @param network.lists The input list
#' @param filter.list A list of metric values to be filtered out. The list 
#' element names have the same names as metric names in network.lists. 
#' The values specify filtering threshold. Everything below the number is 
#' set to 0. Default: NULL
#' @param layout layout to be passed to ggraph. Default: stress. Some useful
#' alternative is dh.
#' 

plot_multi_network = function(network.lists, tissue, filter.list, layout = "stress") {
    
    tissue.nets = get_tissue_dataset_networks(
        tissue,
        network.lists = all.interactions, 
        filter.list = list(MI = 0.2))
    
    for (type in names(tissue.nets)) {
        if (sum(abs(tissue.nets[[type]]) ) == 0  ) {
            tissue.nets[[type]] = NULL
        }
    }
    
    ts.multi.graph = tissue_multi_graph(tissue.nets)
    pp = print_multi_graphs(ts.multi.graph, layout = layout)
    
}


#' Summarize positive and negative interactions across tissues 
#'
#' For each signature pair the function counts how many positive and negative 
#' interactions were observed across tissues and returns a matrix of lists with 
#' two elements called pos and neg or NULL if no interactions were found. 
#' 
#' @param summary.list A list of interactions for individual tissues. Each 
#' element is a matrix.
#' @param tissue.names 
#' @return a list of two matrices for each dataset. The elements of the matrices
#' are lists as elements with two elements called pos and neg for positive and 
#' negative interactions.

summarize_by_tissues = function(summary.list, tissue.names = FALSE) {
    
    all.signatures = sapply(networks.concat, rownames) %>% 
        unlist() %>% 
        unique()
    
    all.mat = matrix(list(), ncol = length(all.signatures),
                                  nrow = length(all.signatures),
                                  dimnames = list(all.signatures, all.signatures))

    for(tissue in names(summary.list)) {
        cat("tissue =", tissue, "\n")
        dt_res = summary.list[[tissue]]
        
        signames = colnames(dt_res)
        
        if (is.null(dt_res)) {
            dt.inds.pos = c()
            dt_inds_neg = c()
            next
        }
        
        dt.inds.pos = which(dt_res > 0, arr.ind = TRUE)
        
        if ( length(dt.inds.pos) != 0) {
            for(i in 1:nrow(dt.inds.pos)) {
                sig1 = signames[dt.inds.pos[i, "row"]]
                sig2 = signames[dt.inds.pos[i, "col"]]
                if (tissue.names) {
                    if (is.null(all.mat[sig1, sig2][[1]]))  {
                        all.mat[sig1, sig2][[1]] = list(pos = "", neg = "")
                    }
                    all.mat[sig1, sig2][[1]][["pos"]] = 
                        paste(tissue, all.mat[sig1, sig2][[1]][["pos"]], sep = ",")
                    
                } else {
                    if (is.null(all.mat[sig1, sig2][[1]]))  {
                        all.mat[sig1, sig2][[1]] = list(pos = 0, neg = 0)
                    }
                    all.mat[sig1, sig2][[1]][["pos"]] = all.mat[sig1, sig2][[1]][["pos"]]  + 1
                }
            }
        }
        
        dt.inds.neg = which(dt_res < 0, arr.ind = TRUE)
        if(length(dt.inds.neg) != 0 ) {
            for(i in 1:nrow(dt.inds.neg)) {
                sig1 = signames[dt.inds.neg[i, "row"]]
                sig2 = signames[dt.inds.neg[i, "col"]]
                if (tissue.names) {
                    if (is.null(all.mat[sig1, sig2][[1]])){
                        all.mat[sig1, sig2][[1]] = list(pos = "", neg = "")
                    }
                    all.mat[sig1, sig2][[1]][["neg"]] = 
                        paste(tissue, all.mat[sig1, sig2][[1]][["neg"]], sep = ",")     
                } else {
                    if (is.null(all.mat[sig1, sig2][[1]])){
                        all.mat[sig1, sig2][[1]] = list(pos = 0, neg = 0)
                    }
                    all.mat[sig1, sig2][[1]][["neg"]] = all.mat[sig1, sig2][[1]][["neg"]] - 1
                }
            }
        }
    }

    return(all.mat)
}
