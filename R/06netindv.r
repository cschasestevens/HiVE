#' Lipid Network Individual Data Plot
#'
#' Extension of hive_base() for plotting individual abundance/expression
#' results onto a base HiVE network. Node/edge colors represent the
#' fold change for the specified comparison for individual lipid
#' species/genes. See HiVE example data for more details.
#' Note that these comparisons are most effective when plotted on
#' a subnetwork.
#'
#' @param type Network type to plot; current lipid networks available are
#' "base" and "oxylipins" (Default), which is a more detailed subnetwork
#' of HiVE base.
#' @param stat_obj A data frame of statistical results containing
#' fold changes and p-values or a Seurat object. The labels must match
#' nodes from the HiVE base network nodes to properly function. Example
#' data are included for formatting help. Can also plot differential
#' gene expression or general expression data if desired.
#' @param col_fc Fold change column name.
#' @param col_p p value column name.
#' @param col_lab Label column for mapping input data onto HiVE nodes.
#' @param col_ct Cell type column name (only applicable if plot_mode is
#' "transcriptomics").
#' @param plot_mode Plot mode; select "lipidomics" for plotting statistical
#' results for lipidomics data or "transcriptomics" for plotting gene expression
#' or DGEA results.
#' @param gex_type Gene expression data type; either "GEX" or "DGEA". Can toggle
#' gex_filt to TRUE to plot only the primary gene between two nodes or FALSE to
#' plot expression/fold changes for all relevant genes between two nodes.
#' @param gex_filt Plot only the primary gene between two nodes?
#' @param tx_size_major Major node label text size.
#' @param tx_size_minor Minor node label text size.
#' @param show_minor_lab Show minor node labels?
#' @param set_incr Set the step distance between expression plot points. The
#' distance is automatically scaled to the number of points if NULL.
#'
#' @return A ggplot2 object plotting individual annotation-level
#' data onto the HiVE base network.
#' @examples
#'
#' # ## Annotation count network
#' # hive_indv(enr_obj = d)
#'
#' @export
hive_indv <- function( # nolint
  type = "oxylipins",
  stat_obj,
  col_fc = "Log2FC",
  col_p = "p.adj",
  col_lab = "label",
  col_ct = "CellType",
  plot_mode = "lipidomics",
  gex_type = "GEX",
  gex_filt = TRUE,
  tx_size_major = 4,
  tx_size_minor = 3,
  show_minor_lab = TRUE,
  set_incr = NULL
) {
  #---- Data input ----
  print("1. Loading input data...")
  if (type == "base") {
    # format edges
    ne <- nedge # nolint
    # format nodes
    nn <- nnode # nolint
    ## gene/enzyme column
    ne[["col.gene"]] <- ne[["enzyme.id"]]
    ## pathway column
    nn[["pathway"]] <- nn[["synthesis.pathway"]]
  }
  if (type == "oxylipins") {
    # format edges
    ne <- edgeoxy # nolint
    # format nodes
    nn <- nodeoxy # nolint
    ## gene/enzyme column
    ne[["col.gene"]] <- ne[["id.gene"]]
    ## pathway column
    nn[["pathway"]] <- nn[["fatty.acid"]]
  }
  # Load stats object
  if (plot_mode == "lipidomics") {
    net1 <- stat_obj
  }
  # Load Seurat if plotting expression level data
  if (plot_mode == "transcriptomics") {
    if (gex_type == "GEX") {
      # Load Seurat and extract all lipid-related genes
      # for the specified set
      net1 <- stat_obj
      net_genes <- unique(unlist(
        lapply(
          seq.int(1, nrow(ne), 1),
          function(i) {
            as.character(
              unlist(strsplit(ne[i, ][["hgnc.symbols"]], ","))
            )
          }
        )
      ))
      net1 <- SeuratObject::FetchData(
        net1,
        vars = c(
          col_ct,
          net_genes[net_genes %in% rownames(net1)]
        )
      )
      if (col_ct != "CellType") {
        net1[["CellType"]] <- net1[[col_ct]]
        net1 <- dplyr::select(
          net1,
          -c(col_ct),
          "CellType", dplyr::everything()
        )
      }
      # calculate mean scaled expression per cell type
      net1[[1]] <- factor(
        as.character(net1[[1]]),
        levels = gtools::mixedsort(unique(as.character(net1[[1]])))
      )
      net1 <- as.data.frame(scale(
        as.matrix(
          magrittr::set_rownames(
            setNames(
              as.data.frame(
                lapply(
                  net1[, 2:ncol(
                    net1
                  )],
                  function(x) {
                    dplyr::select(
                      aggregate(
                        x,
                        list(
                          net1[, 1]
                        ),
                        FUN = mean
                      ),
                      c(2)
                    )
                  }
                )
              ),
              names(net1[, 2:ncol(net1)])
            ),
            levels(net1[, 1])
          )
        ),
        center = TRUE
      ))
    }
    if (gex_type == "DGEA") {
      # Load Seurat and extract all lipid-related genes
      # for the specified set
      net1 <- stat_obj
      net_genes <- unique(unlist(
        lapply(
          seq.int(1, nrow(ne), 1),
          function(i) {
            as.character(
              unlist(strsplit(ne[i, ][["hgnc.symbols"]], ","))
            )
          }
        )
      ))
      net1 <- net1[
        net1[["GENE"]] %in% net_genes,
        c(col_ct, "GENE", "log2FC")
      ]
      # Convert cell type to factor and format as data frame
      net1[[1]] <- factor(
        as.character(net1[[1]]),
        levels = gtools::mixedsort(unique(as.character(net1[[1]])))
      )
      if (col_ct != "CellType") {
        net1[["CellType"]] <- net1[[col_ct]]
        net1 <- dplyr::select(net1, -c(col_ct))
      }
      net1 <- reshape2::dcast(
        net1, CellType ~ GENE, value.var = "log2FC"
      )
      net1[is.na(net1)] <- 0
      net1 <- dplyr::select(
        magrittr::set_rownames(
          net1, levels(net1[["CellType"]])
        ),
        -c("CellType")
      )
    }
  }
  print("2. Setting node attributes...")
  #---- Node attributes ----
  if (plot_mode == "lipidomics") {
    # For each node, set attributes for each individual
    # lipid species
    np_nodes <- setNames(lapply(
      seq.int(1, length(nn[["Node"]]), 1),
      function(i) {
        d_out <- tryCatch(
          {
            # return matches for each node from input stats
            np1 <- net1[net1[[col_lab]] == nn[["Node"]][[i]], ]
            # Check synonyms column to verify that no matches exist
            if (nrow(np1) == 0) {
              np1 <- net1[net1[[col_lab]] == nn[["Synonyms"]][[i]], ]
            }
            if (nrow(np1) == 0) {
              np1 <- paste(
                "No matches for ", nn[["Node"]][[i]],
                " in input data...",
                sep = ""
              )
            }
            # Set node attributes
            if (nrow(np1) > 0) {
              # Set x-y coordinates relative to node position
              ## define function for plotting points around central node
              fun_coord <- function(h, k, r, n) {
                ang1 <- (2 * pi) / n
                d_out <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, n, 1),
                    function(i) {
                      theta <- i * ang1
                      data.frame(
                        "x" = h + r * cos(theta),
                        "y" = k + r * sin(theta)
                      )
                    }
                  )
                )
                return(d_out) # nolint
              }
              if (nrow(np1) == 1) {
                np1 <- cbind(
                  np1,
                  data.frame(
                    "x" = c(nn[i, ][["hive_x"]]),
                    "y" = c(nn[i, ][["hive_y"]])
                  )
                )
              }
              if (nrow(np1) == 2) {
                np1 <- cbind(
                  np1,
                  data.frame(
                    "x" = c(
                      nn[i, ][["hive_x"]] - 0.25,
                      nn[i, ][["hive_x"]] + 0.25
                    ),
                    "y" = c(
                      nn[i, ][["hive_y"]],
                      nn[i, ][["hive_y"]]
                    )
                  )
                )
              }
              if (nrow(np1) >= 3) {
                np1 <- cbind(
                  np1,
                  fun_coord(
                    nn[i, ][["hive_x"]],
                    nn[i, ][["hive_y"]],
                    0.3,
                    nrow(np1)
                  )
                )
              }
              # Set node size and color to emphasize
              # significant points
              np1[["size"]] <- ifelse(np1[[col_p]] < 0.05, 2, 1.5)
              np1[["col.node"]] <- ifelse(np1[[col_p]] < 0.05, "a", "b")
              return(np1) # nolint
            }
          },
          error = function(e) {
            print(
              paste(
                "No matches for ", nn[["Node"]][[i]],
                " in input data...",
                sep = ""
              )
            )
          }
        )
        return(d_out) # nolint
      }
    ), nn[["Node"]])
    np_nodes <- np_nodes[lengths(np_nodes) > 1]
    np_nodes <- dplyr::bind_rows(np_nodes)
    # Scale fold change from 0 to 1 for node fill
    np_nodes[["fill.node"]] <- (
      np_nodes[[col_fc]] - min(np_nodes[[col_fc]], na.rm = TRUE)
    ) / (
      max(np_nodes[[col_fc]], na.rm = TRUE) -
        min(np_nodes[[col_fc]], na.rm = TRUE)
    )
    np_nodes[["fill.node"]] <- ifelse(
      np_nodes[[col_fc]] == 0, 0.5,
      np_nodes[["fill.node"]]
    )
    ### define color scale
    col_fun <- colorRampPalette(col_grad(scm = 5))
    col_len <- col_fun(length(np_nodes[["fill.node"]]))
    np_nodes[["fill.node.names"]] <- col_len
  }
  if (plot_mode == "transcriptomics") {
    print("Plot mode is transcriptomics; using HiVE default node parameters.")
  }
  #---- Edge attributes -----
  print("3. Setting edge attributes...")
  if (plot_mode == "lipidomics") {
    print("Plot mode is lipidomics; using HiVE default edge parameters.")
  }
  if (plot_mode == "transcriptomics") {
    # For each edge, set attributes for each individual
    # lipid-related gene
    np_edges <- setNames(lapply(
      seq.int(1, length(ne[["col.gene"]]), 1),
      function(i) {
        d_out <- tryCatch(
          {
            # return matches for each edge from input data
            np_match <- unlist(strsplit(ne[i, ][["hgnc.symbols"]], ","))
            np1 <- net1[
              ,
              names(net1) %in%
                np_match
            ]
            if (length(np_match[np_match %in% names(net1)]) == 1) {
              np1 <- magrittr::set_rownames(
                setNames(
                  as.data.frame(np1),
                  np_match[np_match %in% names(net1)]
                ),
                rownames(net1)
              )
            }
            if (ncol(np1) == 0) {
              np1 <- paste(
                "No matches for ",
                np_match,
                " in input data...",
                sep = ""
              )
            }
            # Set edge attributes
            if (ncol(np1) > 0) {
              # melt from wide to long
              np2 <- reshape2::melt(
                dplyr::mutate(np1, "CellType" = row.names(np1)),
                id.vars = "CellType"
              )
              # Set x-y coordinates relative to edge position
              ## define function for plotting points around central edge
              ## Set all to the midpoint between
              ## node A ("from") and node B ("to")
              if (nrow(np1) == 1) {
                mid1 <- data.frame(
                  "x" = c(
                    nn[nn[["Node"]] %in% ne[i, ][["from"]], ][["hive_x"]],
                    nn[nn[["Node"]] %in% ne[i, ][["to"]], ][["hive_x"]]
                  ),
                  "y" = c(
                    nn[nn[["Node"]] %in% ne[i, ][["from"]], ][["hive_y"]],
                    nn[nn[["Node"]] %in% ne[i, ][["to"]], ][["hive_y"]]
                  )
                )
                mid1 <- data.frame(
                  "x" = (mid1[["x"]][[1]] + mid1[["x"]][[2]]) / 2,
                  "y" = (mid1[["y"]][[1]] + mid1[["y"]][[2]]) / 2
                )
                np1 <- cbind(
                  np1,
                  mid1
                )
              }
              if (nrow(np1) == 2) {
                mid1 <- data.frame(
                  "x" = c(
                    nn[nn[["Node"]] %in% ne[i, ][["from"]], ][["hive_x"]],
                    nn[nn[["Node"]] %in% ne[i, ][["to"]], ][["hive_x"]]
                  ),
                  "y" = c(
                    nn[nn[["Node"]] %in% ne[i, ][["from"]], ][["hive_y"]],
                    nn[nn[["Node"]] %in% ne[i, ][["to"]], ][["hive_y"]]
                  )
                )
                mid1 <- data.frame(
                  "x" = (mid1[["x"]][[1]] + mid1[["x"]][[2]]) / 2,
                  "y" = (mid1[["y"]][[1]] + mid1[["y"]][[2]]) / 2
                )
                np1 <- cbind(
                  np1,
                  data.frame(
                    "x" = c(
                      mid1[["x"]] - mid1[["x"]] / 2,
                      mid1[["x"]] + mid1[["x"]] / 2
                    ),
                    "y" = c(
                      mid1[["y"]],
                      mid1[["y"]]
                    )
                  )
                )
              }
              if (nrow(np1) >= 3) {
                mid1 <- data.frame(
                  "x" = c(
                    nn[nn[["Node"]] %in% ne[i, ][["from"]], ][["hive_x"]],
                    nn[nn[["Node"]] %in% ne[i, ][["to"]], ][["hive_x"]]
                  ),
                  "y" = c(
                    nn[nn[["Node"]] %in% ne[i, ][["from"]], ][["hive_y"]],
                    nn[nn[["Node"]] %in% ne[i, ][["to"]], ][["hive_y"]]
                  )
                )
                mid1 <- data.frame(
                  "x" = (mid1[["x"]][[1]] + mid1[["x"]][[2]]) / 2,
                  "y" = (mid1[["y"]][[1]] + mid1[["y"]][[2]]) / 2
                )
                # create individual points for all observations around midpoint
                ## set x and y coord increments
                if (is.null(set_incr)) {
                  if (nrow(np1) > 10) {
                    num_x <- (1 / nrow(np1)) / 2
                    num_y <- (1 / nrow(np1)) / 2
                  }
                  if (nrow(np1) < 10 && nrow(np1) > 5) {
                    num_x <- (1 / nrow(np1)) / 4
                    num_y <- (1 / nrow(np1)) / 4
                  }
                  if (nrow(np1) < 5) {
                    num_x <- (1 / nrow(np1)) / 4
                    num_y <- (1 / nrow(np1)) / 4
                  }
                }
                if (!is.null(set_incr)) {
                  num_x <- set_incr
                  num_y <- set_incr
                }
                ## calculate coords for all points
                ### if rows and columns are even in length
                if (nrow(np1) %% 2 == 0 && ncol(np1) %% 2 == 0) {
                  d_out <- data.frame(
                    # determined by cell type
                    "x" = rep(
                      seq(
                        from = mid1[["x"]] -
                          (floor(nrow(np1) / 2) - 1) * num_x,
                        to = mid1[["x"]] +
                          (floor(nrow(np1) / 2)) * num_x,
                        by = num_x
                      ),
                      ncol(np1)
                    ),
                    # determined by gene
                    ## fix error with even number of row/column numbers
                    "y" = sort(rep(
                      seq(
                        from = mid1[["y"]] -
                          (floor(ncol(np1) / 2) - 1) * num_y,
                        to = mid1[["y"]] +
                          (floor(ncol(np1) / 2)) * num_y,
                        by = num_y
                      ),
                      nrow(np1)
                    ))
                  )
                }
                ### if only rows are even
                if (nrow(np1) %% 2 == 0 && ncol(np1) %% 2 == 1) {
                  d_out <- data.frame(
                    # determined by cell type
                    "x" = rep(
                      seq(
                        from = mid1[["x"]] -
                          (floor(nrow(np1) / 2) - 1) * num_x,
                        to = mid1[["x"]] +
                          (floor(nrow(np1) / 2)) * num_x,
                        by = num_x
                      ),
                      ncol(np1)
                    ),
                    # determined by gene
                    ## fix error with even number of row/column numbers
                    "y" = sort(rep(
                      seq(
                        from = mid1[["y"]] -
                          (floor(ncol(np1) / 2)) * num_y,
                        to = mid1[["y"]] +
                          (floor(ncol(np1) / 2)) * num_y,
                        by = num_y
                      ),
                      nrow(np1)
                    ))
                  )
                }
                ### if only columns are even
                if (nrow(np1) %% 2 == 1 && ncol(np1) %% 2 == 0) {
                  d_out <- data.frame(
                    # determined by cell type
                    "x" = rep(
                      seq(
                        from = mid1[["x"]] -
                          (floor(nrow(np1) / 2)) * num_x,
                        to = mid1[["x"]] +
                          (floor(nrow(np1) / 2)) * num_x,
                        by = num_x
                      ),
                      ncol(np1)
                    ),
                    # determined by gene
                    ## fix error with even number of row/column numbers
                    "y" = sort(rep(
                      seq(
                        from = mid1[["y"]] -
                          (floor(ncol(np1) / 2) - 1) * num_y,
                        to = mid1[["y"]] +
                          (floor(ncol(np1) / 2)) * num_y,
                        by = num_y
                      ),
                      nrow(np1)
                    ))
                  )
                }
                ### if neither rows nor columns are even
                if (nrow(np1) %% 2 == 1 && ncol(np1) %% 2 == 1) {
                  d_out <- data.frame(
                    # determined by cell type
                    "x" = rep(
                      seq(
                        from = mid1[["x"]] -
                          (floor(nrow(np1) / 2)) * num_x,
                        to = mid1[["x"]] +
                          (floor(nrow(np1) / 2)) * num_x,
                        by = num_x
                      ),
                      ncol(np1)
                    ),
                    # determined by gene
                    ## fix error with even number of row/column numbers
                    "y" = sort(rep(
                      seq(
                        from = mid1[["y"]] -
                          (floor(ncol(np1) / 2)) * num_y,
                        to = mid1[["y"]] +
                          (floor(ncol(np1) / 2)) * num_y,
                        by = num_y
                      ),
                      nrow(np1)
                    ))
                  )
                }
                np1 <- cbind(
                  np2,
                  d_out
                )
                # Set edge ID column
                np1[["ID"]] <- paste(
                  ne[i, ][["from"]],
                  ne[i, ][["to"]],
                  sep = "_"
                )
              }
              return(np1) # nolint
            }
          },
          error = function(e) {
            print(
              paste(
                "No matches for ",
                np_match,
                " in input data...",
                sep = ""
              )
            )
          }
        )
        return(d_out) # nolint
      }
    ), paste(ne[["from"]], ne[["to"]], sep = "_"))
    np_select <- unlist(lapply(
      seq.int(1, length(np_edges), 1),
      function(j) {
        class(np_edges[[j]]) == "data.frame"
      }
    ))
    np_edges <- np_edges[np_select]
    np_edges <- dplyr::bind_rows(np_edges)
    # Set row and column label coordinates
    ## set column label coords
    d_out2 <- dplyr::bind_rows(
      lapply(
        seq.int(1, length(unique(np_edges[["ID"]])), 1),
        function(i) {
          dout2a <- np_edges[
            np_edges[["ID"]] == unique(np_edges[["ID"]])[[i]],
          ]
          if (length(unique(dout2a[["variable"]])) == 1) {
            dout2b <- data.frame(
              "lab_col_x" = unique(dout2a[["x"]]),
              "lab_col_y" = rep(
                min(dout2a[["y"]]) -
                  ((1 / length(unique(dout2a[["y"]]))) / 10),
                length(unique(dout2a[["x"]]))
              ),
              "CellType" = unique(dout2a[["CellType"]]),
              "ID" = rep(
                unique(dout2a[["ID"]]),
                length(unique(dout2a[["x"]]))
              )
            )
          }
          if (length(unique(dout2a[["variable"]])) == 2) {
            dout2b <- data.frame(
              "lab_col_x" = unique(dout2a[["x"]]),
              "lab_col_y" = rep(
                min(dout2a[["y"]]) -
                  ((1 / length(unique(dout2a[["y"]]))) / 4),
                length(unique(dout2a[["x"]]))
              ),
              "CellType" = unique(dout2a[["CellType"]]),
              "ID" = rep(
                unique(dout2a[["ID"]]),
                length(unique(dout2a[["x"]]))
              )
            )
          }
          if (length(unique(dout2a[["variable"]])) > 2) {
            dout2b <- data.frame(
              "lab_col_x" = unique(dout2a[["x"]]),
              "lab_col_y" = rep(
                min(dout2a[["y"]]) -
                  ((1 / length(unique(dout2a[["y"]]))) / 2),
                length(unique(dout2a[["x"]]))
              ),
              "CellType" = unique(dout2a[["CellType"]]),
              "ID" = rep(
                unique(dout2a[["ID"]]),
                length(unique(dout2a[["x"]]))
              )
            )
          }
          return(dout2b) # nolint
        }
      )
    )
    # row name coords
    d_out3 <- dplyr::bind_rows(
      lapply(
        seq.int(1, length(unique(np_edges[["ID"]])), 1),
        function(i) {
          dout2a <- np_edges[
            np_edges[["ID"]] == unique(np_edges[["ID"]])[[i]],
          ]
          dout2b <- data.frame(
            "lab_row_x" = rep(
              min(dout2a[["x"]]) -
                ((1 / length(unique(dout2a[["x"]]))) / 10),
              length(unique(dout2a[["y"]]))
            ),
            "lab_row_y" = unique(dout2a[["y"]]),
            "GENE" = unique(dout2a[["variable"]]),
            "ID" = rep(
              unique(dout2a[["ID"]]),
              length(unique(dout2a[["y"]]))
            )
          )
          return(dout2b) # nolint
        }
      )
    )
    # Scale fold change from 0 to 1 for node fill
    np_edges[["fill.edge"]] <- (
      np_edges[["value"]] - min(np_edges[["value"]], na.rm = TRUE)
    ) / (
      max(np_edges[["value"]], na.rm = TRUE) -
        min(np_edges[["value"]], na.rm = TRUE)
    )
    np_edges[["fill.edge"]] <- ifelse(
      np_edges[["value"]] == 0, 0.5,
      np_edges[["fill.edge"]]
    )
    ### define color scale
    col_fun <- colorRampPalette(col_grad(scm = 5))
    col_len <- col_fun(length(np_edges[["fill.edge"]]))
    np_edges[["fill.edge.names"]] <- col_len
    ### Set boundary box sizes
    np_edge_box <- lapply(
      unique(np_edges[["ID"]]),
      function(i) {
        d1 <- np_edges[np_edges[["ID"]] == i, ]
        d2 <- dplyr::bind_rows(
          lapply(
            seq.int(1, nrow(d1), 1),
            function(j) {
              d1a <- d1[j, c("x", "y")]
              d2a <- data.frame(
                "x" = c(
                  d1a[["x"]] - 0.025,
                  d1a[["x"]] + 0.025,
                  d1a[["x"]] + 0.025,
                  d1a[["x"]] - 0.025
                ),
                "y" = c(
                  d1a[["y"]] + 0.025,
                  d1a[["y"]] + 0.025,
                  d1a[["y"]] - 0.025,
                  d1a[["y"]] - 0.025
                )
              )
              d2a[["pw"]] <- i
              return(d2a) # nolint
            }
          )
        )
        d3 <- sf::st_as_sf(d2, coords = c("x", "y"))
        d4 <- concaveman::concaveman(d3, concavity = 1)
        d4[["pw"]] <- i
        return(d4) # nolint
      }
    )
    np_edge_box <- do.call(rbind, np_edge_box)
  }
  #---- Create HiVE network ----
  print("4. Generating base network...")
  np <- igraph::graph_from_data_frame(
    ne,
    vertices = nn
  )
  ## Network layout
  npl <- ggraph::create_layout(np, layout = "stress", circular = FALSE)
  # Set manual coordinates for each node
  npl[["x"]] <- npl[["hive_x"]]
  npl[["y"]] <- npl[["hive_y"]]
  if (plot_mode == "lipidomics") {
    npl[["r"]] <- ifelse(
      npl[["name"]] %in%
        unique(np_nodes[[col_lab]]) |
        npl[["Synonyms"]] %in%
          unique(np_nodes[[col_lab]]),
      0.4, 0
    )
  }
  # Define boundaries for each pathway and create node map
  npw <- lapply(
    unique(npl[["pathway"]]),
    function(i) {
      d1 <- npl[npl[["pathway"]] == i, ]
      d2 <- dplyr::bind_rows(
        lapply(
          seq.int(1, nrow(d1), 1),
          function(j) {
            d1a <- d1[j, c("x", "y")]
            d2a <- data.frame(
              "x" = c(
                d1a[["x"]] - 0.5,
                d1a[["x"]] + 0.5,
                d1a[["x"]] + 0.5,
                d1a[["x"]] - 0.5
              ),
              "y" = c(
                d1a[["y"]] + 0.5,
                d1a[["y"]] + 0.5,
                d1a[["y"]] - 0.5,
                d1a[["y"]] - 0.5
              )
            )
            d2a[["pw"]] <- i
            return(d2a) # nolint
          }
        )
      )
      d3 <- sf::st_as_sf(d2, coords = c("x", "y"))
      d4 <- concaveman::concaveman(d3, concavity = 1)
      d4[["pw"]] <- i
      return(d4) # nolint
    }
  )
  npw <- do.call(rbind, npw)
  #---- Plot network ----
  print("5. Plotting HiVE base network")
  if (plot_mode == "lipidomics") {
    # Plot base network
    np_base <- ggraph::ggraph(npl) +
      # color scheme
      ggplot2::scale_color_manual(
        "Pathway",
        values = col_univ() # nolint
      ) +
      ggplot2::scale_fill_manual(
        "Pathway",
        values = col_univ()
      ) +
      # Subclass shading
      ggplot2::geom_sf(
        data = npw,
        ggplot2::aes(
          fill = .data[["pw"]]
        ),
        color = "grey25",
        alpha = 1
      ) +
      # graph edges and attributes
      ggraph::geom_edge_arc(
        ggplot2::aes(label = .data[["col.gene"]]),
        strength = 0.1,
        label_dodge = ggplot2::unit(2, "mm"),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(4, "mm"), type = "closed"
        ),
        start_cap = ggraph::circle(3, "mm"),
        end_cap = ggraph::circle(3, "mm"),
        angle_calc = "along",
        alpha = 0.5,
        label_alpha = 0.5,
        color = "grey50"
      ) +
      ggnewscale::new_scale_fill() +
      ggforce::geom_circle(
        ggplot2::aes(
          x0 = .data[["x"]],
          y0 = .data[["y"]],
          r = .data[["r"]]
        ),
        color = "grey25",
        fill = "grey",
        show.legend = FALSE,
        alpha = 0.9
      )
  }
  if (plot_mode == "transcriptomics") {
    # Plot base network
    np_base <- ggraph::ggraph(npl) +
      # color scheme
      ggplot2::scale_color_manual(
        "Pathway",
        values = col_univ() # nolint
      ) +
      ggplot2::scale_fill_manual(
        "Pathway",
        values = col_univ()
      ) +
      # Subclass shading
      ggplot2::geom_sf(
        data = npw,
        ggplot2::aes(
          fill = .data[["pw"]]
        ),
        color = "grey25",
        alpha = 1
      ) +
      # graph edges and attributes
      ggraph::geom_edge_arc(
        ggplot2::aes(label = .data[["col.gene"]]),
        strength = 0.1,
        label_dodge = ggplot2::unit(2, "mm"),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(4, "mm"), type = "closed"
        ),
        start_cap = ggraph::circle(3, "mm"),
        end_cap = ggraph::circle(3, "mm"),
        angle_calc = "along",
        alpha = 0.5,
        label_alpha = 0.5,
        color = "grey50"
      ) +
      ggnewscale::new_scale_fill() +
      # Edge expression point boundary boxes
      ggplot2::geom_sf(
        data = np_edge_box,
        ggplot2::aes(
          fill = .data[["pw"]]
        ),
        color = "grey25",
        fill = "grey75",
        alpha = 0.9
      )
  }
  #---- Add node attributes ----
  print("6. Plotting node/edge attributes onto network...")
  if (plot_mode == "lipidomics") {
    # Add nodes based on saturation
    np_enr_nodes <- np_base +
      ggnewscale::new_scale_fill() +
      # node points
      ggraph::geom_node_point(
        data = np_nodes,
        ggplot2::aes(
          x = np_nodes[["x"]],
          y = np_nodes[["y"]],
          fill = np_nodes[["fill.node"]]
        ),
        shape = 21,
        size = 5,
        color = "black",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_continuous(
        palette = col_grad(scm = 5)
      )
  }
  if (plot_mode == "transcriptomics") {
    # Add node points
    np_enr_nodes <- np_base +
      # node points
      ggraph::geom_node_point(
        shape = 21,
        size = 5,
        color = "black",
        fill = "grey25",
        show.legend = FALSE
      ) +
      ggnewscale::new_scale_fill() +
      # edge points
      ggraph::geom_node_point(
        data = np_edges,
        ggplot2::aes(
          x = np_edges[["x"]],
          y = np_edges[["y"]],
          fill = np_edges[["fill.edge"]]
        ),
        shape = 22,
        size = 2,
        color = "black",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_continuous(
        palette = col_grad(scm = 5)
      )
  }
  #---- Add all node labels ----
  print("7. Plotting node labels...")
  if (plot_mode == "lipidomics") {
    # Add node labels
    if (show_minor_lab == TRUE) {
      np_enr_labs <- np_enr_nodes +
        ggnewscale::new_scale_color() +
        # Primary node labels
        ggplot2::geom_label(
          ggplot2::aes(
            x = .data[["x"]],
            y = .data[["y"]],
            label = .data[["name"]],
            alpha = ifelse(
              .data[["name"]] %in%
                unique(np_nodes[[col_lab]]) |
                .data[["Synonyms"]] %in%
                  unique(np_nodes[[col_lab]]),
              1, 0.75
            ),
            linewidth = ifelse(
              .data[["name"]] %in%
                unique(np_nodes[[col_lab]]) |
                .data[["Synonyms"]] %in%
                  unique(np_nodes[[col_lab]]),
              0.15, NA
            ),
            nudge_y = ifelse(
              .data[["name"]] %in%
                unique(np_nodes[[col_lab]]) |
                .data[["Synonyms"]] %in%
                  unique(np_nodes[[col_lab]]),
              0.4, 0
            )
          ),
          size = tx_size_major,
          show.legend = FALSE
        ) +
        # Individual node labels
        shadowtext::geom_shadowtext(
          data = np_nodes,
          ggplot2::aes(
            x = np_nodes[["x"]],
            y = np_nodes[["y"]],
            label = np_nodes[["Name"]]
          ),
          nudge_y = 0.1,
          size = tx_size_minor,
          color = "white",
          bg.color = "grey25",
          show.legend = FALSE
        )
    }
    if (show_minor_lab == FALSE) {
      np_enr_labs <- np_enr_nodes +
        ggnewscale::new_scale_color() +
        # Primary node labels
        ggplot2::geom_label(
          ggplot2::aes(
            x = .data[["x"]],
            y = .data[["y"]],
            label = .data[["name"]],
            alpha = ifelse(
              .data[["name"]] %in%
                unique(np_nodes[[col_lab]]) |
                .data[["Synonyms"]] %in%
                  unique(np_nodes[[col_lab]]),
              1, 0.75
            ),
            linewidth = ifelse(
              .data[["name"]] %in%
                unique(np_nodes[[col_lab]]) |
                .data[["Synonyms"]] %in%
                  unique(np_nodes[[col_lab]]),
              0.15, NA
            ),
            nudge_y = ifelse(
              .data[["name"]] %in%
                unique(np_nodes[[col_lab]]) |
                .data[["Synonyms"]] %in%
                  unique(np_nodes[[col_lab]]),
              0.4, 0
            )
          ),
          size = tx_size_major,
          show.legend = FALSE
        )
    }
  }
  if (plot_mode == "transcriptomics") {
    np_enr_labs <- np_enr_nodes +
      # Primary node labels
      ggplot2::geom_label(
        ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          label = .data[["name"]]
        ),
        size = 4,
        color = "white",
        fill = "grey10",
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = d_out3,
        ggplot2::aes(
          x = d_out3[["lab_row_x"]],
          y = d_out3[["lab_row_y"]],
          label = d_out3[["GENE"]]
        ),
        size = 2,
        hjust = 1,
        color = "grey10"
      ) +
      ggplot2::geom_text(
        data = d_out2,
        ggplot2::aes(
          x = d_out2[["lab_col_x"]],
          y = d_out2[["lab_col_y"]],
          label = d_out2[["CellType"]]
        ),
        size = 2,
        angle = 90,
        hjust = 1,
        color = "grey10"
      )
  }
  #---- Add HiVE theme ----
  print("8. Adding HiVE theme...")
  if (
    plot_mode == "lipidomics" ||
      plot_mode == "transcriptomics"
  ) {
    np_enr_thm <- np_enr_labs +
      ggplot2::scale_linewidth(
        range = c(0, 0.15)
      ) +
      # plot theme
      gen_theme() + # nolint
      net_theme() # nolint
  }
  print("HiVE network successfully created!")
  return(np_enr_thm)
}
