#' Lipid Network Enrichment Plot
#'
#' Extension of hive_base() for plotting enrichment results onto
#' a base HiVE network. Node sizes are proportional to the number
#' of individual lipids present in the dataset for each class, node
#' colors represent the proportion of increased/decreased lipids
#' within a specific class, and yellow boxes indicate classes that
#' are significantly enriched in the selected comparison. See
#' HiVE example data for more details.
#'
#' @param level Base level to plot (either 0 or 1). Level 1 contains
#' finer class-based separation while 0 (default) plots higher level
#' lipid classes.
#' @param enr_obj An enrichment object including lipid class labels
#' and the number of lipids within each class. The labels must match
#' nodes from the HiVE base network nodes to properly function. Example
#' data are included for formatting help.
#' @param comp_sat Compare enrichment results split by lipid saturation.
#' This is currently the only option but more will be supported in
#' furture package versions.
#'
#' @return A ggplot2 object plotting a HiVE enrichment network.
#' @examples
#'
#' # ## Annotation count network
#' # hive_enrch(enr_obj = d)
#'
#' @export
hive_enrch <- function( # nolint
  level = 0,
  enr_obj,
  comp_sat = TRUE
) {
  # format edges
  ne <- nedge[nedge[["Level"]] == level, ] # nolint
  # format nodes
  nn <- nnode[nnode[["Level"]] == level, ] # nolint
  # format data
  net1 <- enr_obj
  #---- Enrichment network ----
  print("1. Mapping enrichment results to HiVE network...")
  # Map dataset annotations
  net_anno <- setNames(aggregate(
    unique(net1[, c("label", "n")])[["n"]], # nolint
    list(unique(net1[, c("label", "n")])[["label"]]),
    function(x) sum(x)
  ), c("Node", "cluster_size"))
  ## Adjust node attributes
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
    # define boundary polygon for each subclass
    npw <- dplyr::bind_rows(
      lapply(
        unique(npl[["synthesis.pathway"]]),
        function(i) {
          d1 <- npl[npl[["synthesis.pathway"]] == i, ]
          y1 <- d1[d1[["x"]] < quantile(d1[["x"]], 0.5), ]
          y1 <- y1[order(y1[["y"]]), ]
          y2 <- d1[d1[["x"]] > quantile(d1[["x"]], 0.5), ]
          y2 <- y2[order(y2[["y"]], decreasing = TRUE), ]
          d4 <- dplyr::bind_rows(y1, y2)[, c("x", "y")]
          d4[["pw"]] <- i
          d2 <- data.frame(
            "x" = c(
              min(d1[["x"]]) - abs(0.15 * min(d1[["x"]])),
              min(d1[["x"]]) - abs(0.15 * min(d1[["x"]])),
              max(d1[["x"]]) + abs(0.15 * max(d1[["x"]])),
              max(d1[["x"]]) + abs(0.15 * max(d1[["x"]]))
            ),
            "y" = c(
              min(d1[["y"]]) - abs(0.15 * min(d1[["y"]])),
              max(d1[["y"]]) + abs(0.15 * max(d1[["y"]])),
              max(d1[["y"]]) + abs(0.15 * max(d1[["y"]])),
              min(d1[["y"]]) - abs(0.15 * min(d1[["y"]]))
            ),
            "pw" = i
          )
          return(d2) # nolint
        }
      )
    )
    npw[["pw"]] <- factor(
      npw[["pw"]],
      levels = c(
        "Phospholipid Biosynthesis", "Fatty Acid Biosynthesis",
        "Lipid Transport", "Plasmalogen Biosynthesis",
        "Oxylipin Biosynthesis", "Sterol Lipid Biosynthesis",
        "Cholesterol Catabolism", "Sphingolipid Biosynthesis",
        "Glycerolipid Biosynthesis"
      )
    )
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
      ggplot2::geom_polygon(
        data = npw,
        ggplot2::aes(
          x = .data[["x"]], # nolint
          y = .data[["y"]],
          fill = .data[["pw"]]
        ),
        color = "grey25",
        alpha = 1
      ) +
      # edge attributes
      ggraph::geom_edge_link(
        ggplot2::aes(
          label = ifelse(
            .data[["set_alpha_edge"]] == 0.25, # nolint
            "",
            .data[["enzyme.id"]]
          ),
          alpha = .data[["set_alpha_edge"]]
        ),
        label_dodge = ggplot2::unit(2, "mm"),
        alpha = 0.5,
        label_alpha = 0.0,
        angle_calc = "along",
        color = "grey50",
        show.legend = FALSE
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
