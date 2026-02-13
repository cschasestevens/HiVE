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
#' @param lab_a Edge label alpha value.
#' @param stat_obj A data frame of statistical results containing
#' fold changes and p-values. The labels must match
#' nodes from the HiVE base network nodes to properly function. Example
#' data are included for formatting help.
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
  type = "base",
  lab_a = 0.5,
  stat_obj
) {
  ## START HERE: add compatibility for Seurat objects
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
  # format data
  net1 <- enr_obj
  #---- Enrichment network ----
  print("1. Mapping enrichment results to HiVE network...")
  # Map dataset annotations
  if (type == "base") {
    net_anno <- stats::setNames(stats::aggregate(
      unique(net1[, c("label", "n")])[["n"]], # nolint
      list(unique(net1[, c("label", "n")])[["label"]]),
      function(x) sum(x)
    ), c("Node", "cluster_size"))
    nn <- dplyr::left_join(
      nn,
      net_anno,
      by = "Node"
    )
    nn[["cluster_size"]] <- ifelse(
      is.na(nn[["cluster_size"]]),
      0,
      nn[["cluster_size"]]
    )
  }
  ## Alternate mapping function if plotting subnetworks
  if (type == "oxylipins") {
    net_anno <- stats::setNames(stats::aggregate(
      net1[, c("label", "n")][["n"]], # nolint
      list(net1[, c("label", "n")][["label"]]),
      function(x) sum(x)
    ), c("Node", "cluster_size"))
    # Join and set cluster_size
    nn <- dplyr::left_join(
      nn,
      net_anno,
      by = "Node"
    )
    nn[["cluster_size"]] <- ifelse(
      is.na(nn[["cluster_size"]]),
      0,
      nn[["cluster_size"]]
    )
    # Check for synonym matches in node list
    nn[nn[["Synonyms"]] %in% net_anno[["Node"]], ][["cluster_size"]] <-
      nn[nn[["Synonyms"]] %in% net_anno[["Node"]], ][["cluster_size"]] +
      net_anno[net_anno[["Node"]] %in% nn[["Synonyms"]], ][["cluster_size"]]
  }
  nn[["set_alpha"]] <- (
    (nn[["cluster_size"]] - min(nn[["cluster_size"]])) /
      (max(nn[["cluster_size"]])) - min(nn[["cluster_size"]])
  )
  nn[["set_clus_size"]] <- (
    (nn[["cluster_size"]] - min(nn[["cluster_size"]])) /
      (max(nn[["cluster_size"]])) - min(nn[["cluster_size"]])
  )
  ## Adjust edge attributes
  ne[["set_alpha_edge"]] <- ifelse(
    ne[["from"]] %in%
      nn[nn[["set_clus_size"]] > 0, "Node"],
    1,
    0.25
  )
  # Map enrichment statistics
  #---- If comparing lipids based on saturation ----
  ## Join enrichment results with network plot
  if (comp_sat == TRUE) {
    net1[["Saturation"]] <- ifelse(
      grepl("\\.Uns|\\.uns", net1[["label.saturation"]]),
      "Unsat",
      ifelse(
        grepl("\\.Sat|\\.sat", net1[["label.saturation"]]),
        "Sat",
        ""
      )
    )
    net_enrch <- as.data.frame(tidyr::pivot_wider(
      net1[, c( # nolint
        "label", "Saturation",
        "Ratio", "FDR", "n"
      )],
      names_from = c(.data[["Saturation"]]), # nolint
      values_from = c(.data[["Ratio"]], .data[["FDR"]], .data[["n"]])
    ))
    net_name1 <- names(net_enrch[, 2:ncol(net_enrch)]) # nolint
    net_enrch[["Node"]] <- net_enrch[["label"]]
    nn <- dplyr::left_join(
      nn,
      net_enrch,
      by = "Node"
    )
    ## Network
    print("2. Creating HiVE network and setting attributes...")
    np <- igraph::graph_from_data_frame(
      ne,
      vertices = nn
    )
    ## Network layout
    npl <- ggraph::create_layout(np, layout = "stress", circular = FALSE)
    # Set manual coordinates for each node
    npl[["x"]] <- npl[["hive_x"]]
    npl[["y"]] <- npl[["hive_y"]]
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
    # Node box color to denote significant clusters
    npl_box_color <- npl[grepl("FDR", names(npl))]
    npl_box_color[is.na(npl_box_color)] <- 1
    npl[["box_color"]] <- unlist(lapply(
      seq.int(1, nrow(npl_box_color), 1),
      function(i) {
        min(npl_box_color[i, ]) < 0.05
      }
    ))
    # Node attributes
    np_nodes <- setNames(
      lapply(
        c("Sat", "Unsat", ""),
        function(j) {
          dnode <- data.frame(
            "x" = npl[["x"]],
            "y" = npl[["y"]],
            "size" = ifelse(
              is.na(npl[[paste(
                "n", j, sep = "_"
              )]]) |
                unlist(
                  lapply(
                    npl[[paste(
                      "n", j, sep = "_"
                    )]],
                    function(x) is.null(x)
                  )
                ),
              0,
              sqrt(as.numeric(unlist(npl[[paste(
                "n", j, sep = "_"
              )]])))
            ),
            "fill" = ifelse(
              is.na(npl[[paste(
                "Ratio", j, sep = "_"
              )]]) |
                unlist(
                  lapply(
                    npl[[paste(
                      "Ratio", j, sep = "_"
                    )]],
                    function(x) is.null(x)
                  )
                ),
              0.5,
              unlist(npl[[paste(
                "Ratio", j, sep = "_"
              )]])
            ),
            "color" = ifelse(
              is.na(npl[[paste(
                "Ratio", j, sep = "_"
              )]]) |
                unlist(
                  lapply(
                    npl[[paste(
                      "Ratio", j, sep = "_"
                    )]],
                    function(x) is.null(x)
                  )
                ),
              "b",
              "a"
            ),
            "lab" = ifelse(
              unlist(
                lapply(
                  npl[[paste("n", j, sep = "_")]],
                  function(x) is.null(x)
                )
              ),
              "",
              j
            ),
            "lab_FDR" = unlist(
              lapply(
                npl[[paste("FDR", j, sep = "_")]],
                function(x) {
                  d1 <- ifelse(
                    !is.null(x),
                    x > 0.05 | is.na(x),
                    TRUE
                  )
                  d1 <- ifelse(
                    d1 == FALSE,
                    paste("FDR = ", round(x, digits = 4)),
                    ""
                  )
                  return(d1) # nolint
                }
              )
            ),
            "lab_size" = ifelse(
              npl[["cluster_size"]] == 0,
              0,
              1.5
            )
          )
          dnode <- dnode[dnode[["size"]] > 0, ]
          return(dnode) # nolint
        }
      ),
      c("Sat", "Unsat", "")
    )
    print("3. Generate base HiVE network...")
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
        ggplot2::aes(
          label = ifelse(
            .data[["set_alpha_edge"]] == 0.25, # nolint
            "",
            .data[["col.gene"]]
          ),
          alpha = .data[["set_alpha_edge"]]
        ),
        curvature = 0.1,
        label_dodge = ggplot2::unit(2, "mm"),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(4, "mm"), type = "closed"
        ),
        start_cap = ggraph::circle(3, "mm"),
        end_cap = ggraph::circle(3, "mm"),
        angle_calc = "along",
        alpha = 0.5,
        label_alpha = lab_a,
        color = "grey50"
      ) +
      ggnewscale::new_scale_fill() +
      # Background for node panels
      ggplot2::geom_rect(
        ggplot2::aes(
          fill = .data[["box_color"]],
          xmin = ifelse(
            .data[["cluster_size"]] == 0,
            0,
            x - 0.35 # nolint
          ),
          ymin = ifelse(
            .data[["cluster_size"]] == 0,
            0,
            y - 0.35 # nolint
          ),
          ymax = ifelse(
            .data[["cluster_size"]] == 0,
            0,
            y + 0.25
          ),
          xmax = ifelse(
            .data[["cluster_size"]] == 0,
            0,
            x + 0.35
          )
        ),
        alpha = 0.9,
        color = "grey25",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        "Significant",
        values = c("grey", "goldenrod1")
      )
    print("4. Add node attributes...")
    # Add nodes based on saturation
    np_enr_nodes <- np_base +
      # node points
      ggraph::geom_node_point(
        ggplot2::aes(
          alpha = ifelse(
            .data[["set_alpha"]] == 0, # nolint
            1,
            0
          )
        ),
        color = "grey75",
        show.legend = FALSE
      ) +
      ## Change scale
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      ## Enrichment nodes
      ## Saturated lipids
      ggraph::geom_node_point(
        data = np_nodes[["Sat"]],
        ggplot2::aes(
          x = np_nodes[["Sat"]][["x"]],
          y = np_nodes[["Sat"]][["y"]],
          size = np_nodes[["Sat"]][["size"]],
          fill = np_nodes[["Sat"]][["fill"]],
          color = np_nodes[["Sat"]][["color"]]
        ),
        position = ggplot2::position_nudge(
          x = -0.25,
          y = -0.1
        ),
        shape = 21,
        show.legend = FALSE
      ) +
      ggraph::geom_node_point(
        data = np_nodes[["Unsat"]],
        ggplot2::aes(
          x = np_nodes[["Unsat"]][["x"]],
          y = np_nodes[["Unsat"]][["y"]],
          size = np_nodes[["Unsat"]][["size"]],
          fill = np_nodes[["Unsat"]][["fill"]],
          color = np_nodes[["Unsat"]][["color"]]
        ),
        position = ggplot2::position_nudge(
          x = 0.25,
          y = -0.1
        ),
        shape = 21,
        show.legend = FALSE
      ) +
      ggraph::geom_node_point(
        data = np_nodes[[3]],
        ggplot2::aes(
          x = np_nodes[[3]][["x"]],
          y = np_nodes[[3]][["y"]],
          size = np_nodes[[3]][["size"]],
          fill = np_nodes[[3]][["fill"]],
          color = np_nodes[[3]][["color"]]
        ),
        position = ggplot2::position_nudge(
          x = 0,
          y = -0.1
        ),
        shape = 21,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(
        "Pathway",
        values = c("black", "white")
      )
    print("5. Add node labels...")
    # Add node labels
    np_enr_labs <- np_enr_nodes +
      ggnewscale::new_scale_color() +
      # ChemRICH node labels
      ## Saturated
      shadowtext::geom_shadowtext(
        data = np_nodes[["Sat"]],
        ggplot2::aes(
          x = np_nodes[["Sat"]][["x"]],
          y = np_nodes[["Sat"]][["y"]],
          label = np_nodes[["Sat"]][["lab"]],
          size = np_nodes[["Sat"]][["lab_size"]]
        ),
        nudge_x = -0.25,
        nudge_y = -0.2,
        color = "white",
        bg.color = "grey25",
        show.legend = FALSE
      ) +
      ## Unsaturated
      shadowtext::geom_shadowtext(
        data = np_nodes[["Unsat"]],
        ggplot2::aes(
          x = np_nodes[["Unsat"]][["x"]],
          y = np_nodes[["Unsat"]][["y"]],
          label = np_nodes[["Unsat"]][["lab"]],
          size = np_nodes[["Unsat"]][["lab_size"]]
        ),
        nudge_x = 0.25,
        nudge_y = -0.2,
        color = "white",
        bg.color = "grey25",
        show.legend = FALSE
      ) +
      ## No saturation
      shadowtext::geom_shadowtext(
        data = np_nodes[[3]],
        ggplot2::aes(
          x = np_nodes[[3]][["x"]],
          y = np_nodes[[3]][["y"]],
          label = np_nodes[[3]][["lab"]],
          size = np_nodes[[3]][["lab_size"]]
        ),
        nudge_x = 0,
        nudge_y = -0.2,
        color = "white",
        bg.color = "grey25",
        show.legend = FALSE
      ) +
      # text labels
      ggrepel::geom_text_repel(
        ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          label = .data[["label"]],
          alpha = ifelse(
            .data[["cluster_size"]] == 0,
            0.1,
            1
          ),
          size = ifelse(
            .data[["cluster_size"]] == 0,
            0.05,
            2
          )
        ),
        color = "white",
        bg.r = 0.03,
        bg.color = "grey10",
        nudge_x = 0,
        nudge_y = 0.1 * 1,
        show.legend = FALSE
      )
    print("6. Add HiVE network theme...")
    # Add theme
    np_enr_thm <- np_enr_labs +
      ggplot2::scale_fill_gradientn(
        "Increased Ratio",
        colors = c(
          "dodgerblue4",
          "azure",
          "firebrick3"
        )
      ) +
      ## Change scale
      ggnewscale::new_scale_color() +
      ggplot2::scale_size_area(
        max_size = 10
      ) +
      ggplot2::scale_color_manual(
        "Pathway",
        values = col_univ()
      ) +
      # plot theme
      gen_theme() + # nolint
      net_theme() # nolint
    print("HiVE enrichment network successfully created!")
  }
  return(np_enr_thm)
}
