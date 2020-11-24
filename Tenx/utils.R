colors_overload <- union(ggthemes::tableau_color_pal('Tableau 20')(20), RColorBrewer::brewer.pal(12, 'Set3'))
colors_overload <- c(colors_overload, 'black')

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}

create_object <- function(
    exprs_raw, meta_data, npcs=30, var_genes=NULL, split_vargenes_by=NULL, nvargenes=2000, verbose=FALSE,
    max_mt=Inf, min_ngene=0, max_numi=Inf, gene_exclude_pattern='^MT-|^RPS|^RPL'
) {
    obj <- list()
    
    message('start filter')
    meta_data <- meta_data %>% 
        dplyr::filter(percent_mito < max_mt & nGene >= min_ngene & nUMI < max_numi)
    exprs_raw <- exprs_raw[, meta_data$CellID]
    
    if (!'weight' %in% colnames(meta_data)) {
        warning('weights not initialized in metadata. Setting all to 1.')
        meta_data$weight <- rep(1, nrow(meta_data))
    }
    obj$meta_data <- meta_data
    obj$exprs_raw <- exprs_raw
    
    message('start normalization')
    t <- system.time({
        genes_use <- which(Matrix::rowSums(exprs_raw != 0) >= 10)
        genes_use <- genes_use[which(!grepl(gene_exclude_pattern, names(genes_use)))]
        exprs_norm <- exprs_raw[genes_use, ] %>% 
#             normalizeData(method = 'log', scaling_factor = median(meta_data$nUMI))        
            normalizeData(method = 'log', scaling_factor = 1e4)   
        obj$exprs_norm <- exprs_norm
    })
    if (verbose) {
        print(t)
        message('Finished normalization')
    }
    
    if (is.null(var_genes)) {
        message('start vargenes')
        t <- system.time({
            if (missing(split_vargenes_by)) {
                var_genes <- singlecellmethods::vargenes_vst(exprs_norm, topn = nvargenes)    
            } else {
                var_genes <- singlecellmethods::vargenes_vst(exprs_norm, meta_data[[split_vargenes_by]], topn = nvargenes) 
            }
        })
        if (verbose) {
            print(t)
            message('Finished vargenes')
        }        
    } else {
        ## for safety
        var_genes <- intersect(var_genes, rownames(obj$exprs_norm))
    }
    
    obj$var_genes <- var_genes

    message('start pca')
    t <- system.time({
        pca_res <- weighted_pca(exprs_norm[obj$var_genes, ], meta_data[['weight']], npc=npcs)
    })
    if (verbose) {
        print(t)
        message('Finished PCA')
    }
    obj$V <- pca_res$embeddings
    obj$loadings <- pca_res$loadings
    obj$vargenes_means_sds <- pca_res$vargenes

#     message('start UMAP')
#     t <- system.time({
#         obj$umap_before_fname <- tempfile(tmpdir = '/data/srlab/ik936/Roche/data/cache', pattern = 'umap_')
#         umap_res <- do_umap(obj$V, obj$umap_before_fname)
#         obj$umap_before <- umap_res$embedding
#         obj$adj_before <- umap_res$adj
#         obj$knn_before <- umap_res$knnadj
#     })    
#     if (verbose) {
#         message('Finished UMAP')
#         print(t)
#     }
    
    return(obj)
}


do_harmony <- function(obj, vars, umap_min_dist=.1, max.iter.cluster = 20, ...) {
    
    ## run Harmony
    hres <- HarmonyMatrix(obj$V, obj$meta_data, vars, 
                          do_pca = FALSE, 
                          max.iter.cluster = max.iter.cluster, 
                          return_object = TRUE, ...)
    
    ## save relevant fields for downstream analysis
    obj$Z_cos <- hres$Z_cos
    obj$Z_corr <- hres$Z_corr
    obj$R <- hres$R
    obj$betas <- harmony:::moe_ridge_get_betas(hres)
    obj$kmeans_rounds <- hres$kmeans_rounds
    obj$objective_kmeans <- hres$objective_kmeans
    obj$use_weights <- hres$use_weights
    obj$weights <- hres$weights
    
    
    ## recompute UMAP on Harmonized PCs
    obj$umap_after_fname <- tempfile(tmpdir = '/data/srlab/ik936/Roche/data/cache', pattern = 'umap_')
    umap_res <- do_umap(t(obj$Z_cos), obj$umap_after_fname)
    obj$umap_after <- umap_res$embedding
    obj$adj_after <- umap_res$adj
    obj$knn_after <- umap_res$knnadj
    return(obj)
}




do_umap <- function(
    Xmat, cache_fname=NULL, 
    .spread=0.3, .min_dist=0.05,
    .metric='euclidean', .init='laplacian',
    .a=NULL, .b=NULL,
    .n_components=2L,
    ...
) {
    umap_object <- uwot::umap(
        X = Xmat,
        n_threads = 20,
        n_neighbors = 30L,
        n_components = .n_components,
        metric = .metric,
        init = .init,
        n_epochs = NULL,
        learning_rate = 1.0,
#         min_dist = 0.3,
#         spread = 1.0,
        min_dist = .min_dist, 
        spread = .spread,        
        set_op_mix_ratio = 1.0,
        local_connectivity = 1L,
        repulsion_strength = 1,
        negative_sample_rate = 1,
        a = .a,
        b = .b,
        fast_sgd = FALSE,
        verbose = FALSE, 
#         ret_model = TRUE,
#         ret_nn = TRUE
        ret_extra = c('nn', 'fgraph', 'model'),
        ...
    ) 

    ## save object for mapping new data
    if (!is.null(cache_fname)) {
        uwot::save_uwot(umap_object, file = cache_fname)#, unload = FALSE, verbose = FALSE)
    }

    ## fxn from dist to kernel from UWOT 
    nn_idx <- umap_object$nn[[1]]$idx
    adj <- Matrix::sparseMatrix(
        i = rep(1:nrow(nn_idx), each = ncol(nn_idx)), 
        j = c(t(nn_idx)), 
        x = c(t(exp(-(pmax(umap_object$nn[[1]]$dist, .min_dist) - .min_dist)/.spread)))
    )
    
    ## return embeddings
    return(list(
        embedding=umap_object$embedding,
        adj=umap_object$fgraph + Matrix::Diagonal(n = nrow(umap_object$fgraph)),
        knnadj=adj
    ))
}


do_cluster <- function(
    obj, adj_name, resolutions,
    force_snn=FALSE,
    append_cols=FALSE,
    do_weights = FALSE, 
    slot_name = 'clusters_df',
    ...
) {
    ## cluster 
    library(singlecellmethods)
#     if (!'snn' %in% names(obj)| force_snn) {        
# #         ifelse(
# #             'Z_cos' %in% names(obj), 
# #             Z_use <- t(obj$Z_cos), 
# #             Z_use <- obj$V
# #         )

#         ## Assumes that KNN already computed (e.g. from UMAP)
#         snn <- Matrix::tcrossprod(obj[[adj_name]])
#         nn_k <- sum(obj[[adj_name]][1, ] > 0)        
#         snn@x <- snn@x / (2 * nn_k - snn@x)
#         obj$snn <- snn %>% as('dgCMatrix') %>% drop0()
        
# #         obj$snn <- singlecellmethods:::buildSNN_fromFeatures(Z_use, prune_snn = 1/25, nn_k = 50, nn_eps = 0)
#         if (do_weights) {
#             obj$snn <- obj$snn %*% Matrix::Diagonal(x = obj$meta_data[['weight']])        
#         }
#     }
#     message('finished SNN')

    ## For now, just always do parallel
    future::plan(multiprocess) 
    
    ## save this separately so as not to pass the full object to future_map
#     snn_use <- obj$snn
    adj_use <- obj[[adj_name]]
    adj_size <- as.integer(pryr::object_size(adj_use))
    if (adj_size > 5e8) {
        options(future.globals.maxSize=1.5*adj_size)
    }
    
    res_new <- future_map(resolutions, function(resolution) {
        message(resolution)
        as.character(Seurat:::RunModularityClustering(adj_use, resolution = resolution, print.output = FALSE, ...))
    }) %>% 
        bind_cols()
    message('finished Louvain')
    
    res_new <- apply(res_new, 2, as.character)    
    if (append_cols) {
        obj[[slot_name]] <- cbind(obj[[slot_name]], res_new)
    } else {
        obj[[slot_name]] <- res_new
    }
    obj[[slot_name]] <- data.frame(obj[[slot_name]])
    
    colnames(obj[[slot_name]]) <- paste0('res', seq(ncol(obj[[slot_name]])))
    
#     ## find markers 
#     obj$markers <- apply(obj[[slot_name]], 2, function(clusters) {
#         wilcoxauc(obj$exprs_norm, clusters)
#     })
#     names(obj$markers) <- paste0('res', seq(length(resolutions)))
    
    return(obj)
}


# name_clusters <- function(obj, cluster_name, new_name, name_list) {
# #     message('TODO: include error checking into name_clusters')
#     clusters <- obj$clusters_df[, cluster_name]
#     cluster_labels <- Reduce(rbind, lapply(names(name_list), function(y) {
#         data.table(cell_type = y, cluster_ids = name_list[[y]])
#     }))
#     res <- data.table(cluster_ids = clusters) %>%
#         dplyr::left_join(cluster_labels, by = "cluster_ids") %>% 
#         dplyr::select(-cluster_ids) %>% 
#         with(cell_type)
    
#     if (length(res) != nrow(obj$meta_data)) {
#         stop('cluster names dont match number of cells in meta_data')
#     }
#     obj$meta_data[new_name] <- res
#     return(obj)
# }

do_scatter <- function (umap_use, meta_data, label_name, facet_var, no_guides = TRUE, 
    do_labels = TRUE, nice_names, palette_use = colors_overload, 
    pt_size = 4, point_size = 0.5, pt_shape = ".", base_size = 20, 
    do_points = TRUE, do_density = FALSE, h = 3, w = 4, 
                        alpha_fore=1, alpha_back=.3, color_back='lightgrey', 
                       nrow = 1, do_raster = FALSE) 
{
    if (do_raster) {
        geom_point_fxn <- function(...) geom_point_rast(..., width = w, height = h)
    } else {
        geom_point_fxn <- geom_point
    }
    
    plt_df <- umap_use %>% data.frame() %>% cbind(meta_data) %>% 
        dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]
    if (!missing(nice_names)) {
        plt_df %<>% dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        plt_df[[label_name]] <- plt_df$nice_name
    }
    
    plt <- plt_df %>% ggplot(aes_string("X1", "X2", col = label_name, 
        fill = label_name)) + 
#         theme_tufte(base_size = base_size) + 
#         theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, 
            alpha = 1, shape = 16, size = 4)), alpha = FALSE) + 
        scale_color_manual(values = palette_use) + scale_fill_manual(values = palette_use) + 
        theme(plot.title = element_text(hjust = 0.5)) + labs(x = "UMAP 1", 
        y = "UMAP 2")
    if (do_points) {
        ## this facets while keeping non-facet points in the background
        if (!missing(facet_var)) {
            if (!is(facet_var, 'quosure')) {
                stop('facet_var must be a quosure. e.g. quo(\'donor\')')
            }            

            plt <- plt + geom_point_fxn(
                data = dplyr::select(plt_df, -!!facet_var), 
                shape = pt_shape, size = point_size,
                color = color_back, fill = color_back, alpha = alpha_back
            ) +
                facet_wrap(vars(!!facet_var), nrow = nrow)
        }
        plt <- plt + geom_point_fxn(shape = pt_shape, size = point_size, alpha = alpha_fore)
    }
    if (do_density) 
        plt <- plt + geom_density_2d()
    if (no_guides) 
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    if (do_labels) {
        plt <- plt + 
#             geom_text_repel(
#                 data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
#                 label.size = NA, aes_string(label = label_name), 
#                 color = "black", 
#                 size = pt_size, alpha = 1, segment.size = 0
#             ) + 
            geom_label(
                data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
                label.size = NA, aes_string(label = label_name, color = label_name), 
#                 color = "black", 
                fill = 'white', 
                size = pt_size, alpha = .6, segment.size = 0
            ) + 
            geom_text(
                data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
                label.size = NA, aes_string(label = label_name, color = label_name), 
#                 color = "black", 
                size = pt_size, alpha = 1, segment.size = 0
            ) + 
            guides(col = FALSE, fill = FALSE)
    }
    return(plt)
}


setupVals <- function(data_mat, feature, qlo, qhi) {
    .x <- data_mat[feature, , drop = FALSE] %>% as("dgTMatrix")
    cutoffs <- quantileSparse(.x, c(qlo, qhi))
    cutoffs[2] <- max(cutoffs[2], min(.x@x))
    if (qlo == 0 & qhi == 1) {
        return(.x)
    } 
    
    if (qlo > 0) {
        .x@x[.x@x < cutoffs[1]] <- cutoffs[1]
#         message(sprintf("For %s, lo = %.3f", feature, ifelse(length(.x@x) == ncol(.x), cutoffs[1], NA)))
    }
    if (qhi < 1) {
        .x@x[.x@x > cutoffs[2]] <- cutoffs[2]
#         message(sprintf("For %s, hi = %.3f", feature, cutoffs[2]))
        
    }
    return(.x)
}


quantileSparse <- function(.x, qlist) {
    ratio_zero <- 1 - (length(.x@x) / ncol(.x))
    q_nz <- which(qlist > ratio_zero)
    q_adj <- (qlist[q_nz] - ratio_zero) / (1 - ratio_zero)
    res <- rep(0, length(qlist))
    res[q_nz] <- quantile(.x@x, q_adj)
    res
}

## TODO: test is feature is present
## TODO: allow for different cutoffs, for each marker
## TODO: somehow draw canvas first, then do plotting? 
library(patchwork)
library(ggthemes)

plotFeatures <- function(data_mat, dim_df, features, nrow = 1, 
                         qlo = 0.05, qhi = 1, order_by_expression = FALSE, 
                         pt_shape = 16, pt_size = .5, no_guide = FALSE,
                         .xlim = c(NA, NA), .ylim = c(NA, NA), color_high = muted("blue")) {
    plt_df <- data.frame(dim_df[, 1:2])
    colnames(plt_df) <- c("X1", "X2")


    plt_list <- lapply(features, function(feature) {
        .x <- setupVals(data_mat, feature, qlo, qhi)
        plt_df$value <- 0
        plt_df[.x@j + 1, "value"] <- .x@x
        if (order_by_expression) {
            plt_df %<>% dplyr::arrange(value)             
        } else {
            plt_df %<>% dplyr::sample_frac(1L)
        }

        plt <- plt_df %>% 
            ggplot(aes(X1, X2, color = value)) + 
#             geom_point_rast(dpi = 300, width = 6, height = 4, size = .5, shape = pt_shape) + 
            geom_point(shape = ".") + 
            scale_color_gradient2(na.value = "lightgrey", mid = "lightgrey", midpoint = 0, high = color_high) + 
#             theme_tufte(base_size = 14, base_family = "Helvetica") + 
#             theme(panel.background = element_rect(), plot.title = element_text(hjust = .5)) +
            theme(plot.title = element_text(hjust = .5)) +
            labs(x = "UMAP 1", y = "UMAP 2", title = feature) + 
            NULL
        if (no_guide) {
            plt <- plt + 
            guides(color = FALSE) 
        }
        
        if (sum(is.na(.xlim)) < 2) 
            plt <- plt + xlim(.xlim)
        if (sum(is.na(.ylim)) < 2) 
            plt <- plt + ylim(.ylim)
        plt

    })
    if (length(plt_list) > 1) {
        Reduce(`+`, plt_list) + patchwork::plot_layout(nrow = nrow)
    } else {
        plt_list[[1]]
    }
}


plot_clusters <- function(obj, umap_use='umap_before', resnames=NULL, slot_name='clusters_df') {
    if (is.null(resnames)) {
        resnames <- colnames(obj[[slot_name]])
    }
    res <- lapply(resnames, function(resname) {
        do_scatter(obj[[umap_use]], obj[[slot_name]], resname, pt_size = 8) + 
            labs(title = resname)
    }) 
    if (length(res) > 1) {
        res <- purrr::reduce(res, `+`) 
    }
    return(res)
}







