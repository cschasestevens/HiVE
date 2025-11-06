#' Lipid Network Annotation Plot
#'
#' Extension of hive_base() for adding annotation count information
#' to a base HiVE network. Node sizes are proportional to the number
#' of individual lipids present in the dataset for each class. See
#' HiVE example data for more details.
#'
#' @param level Base level to plot (either 0 or 1). Level 1 contains
#' finer class-based separation while 0 (default) plots higher level
#' lipid classes.
#' @param enr_obj An enrichment object including lipid class labels
#' and the number of lipids within each class. The labels must match
#' nodes from the HiVE base network nodes to properly function. Example
#' data are included for formatting help.
#'
#' @return A ggplot2 object plotting the HiVE annotation network.
#' @examples
#'
#' # ## Annotation count network
#' # hive_anno(enr_obj = d)
#'
#' @export
hive_anno <- function(
  level = 0,
  enr_obj
) {
  # format edges
  ne <- nedge[nedge[["Level"]] == level, ] # nolint
  # format nodes
  nn <- nnode[nnode[["Level"]] == level, ] # nolint
  # format data
  net1 <- enr_obj
  #---- Annotation network ----
  # Map dataset annotations
  net_anno <- stats::setNames(stats::aggregate(
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
  # Network
  np <- igraph::graph_from_data_frame(
    ne,
    vertices = nn
  )
  npl <- ggraph::create_layout(np, layout = "stress", circular = FALSE)
  # define boundary polygon for each subclass
  npw <- dplyr::bind_rows(
    lapply(
      unique(npl[["synthesis.pathway"]]),
      function(i) {
        d1 <- npl[npl[["synthesis.pathway"]] == i, ]
        y1 <- d1[d1[["x"]] < stats::quantile(d1[["x"]], 0.5), ]
        y1 <- y1[order(y1[["y"]]), ]
        y2 <- d1[d1[["x"]] > stats::quantile(d1[["x"]], 0.5), ]
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
  # Plot Annotations
  ## Plot
  np_anno <- ggraph::ggraph(npl) + # Base graph
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
    # node points
    ggraph::geom_node_point(
      ggplot2::aes(
        alpha = ifelse(
          .data[["set_alpha"]] == 0, # nolint
          1,
          0
        )
      ),
      color = "grey50",
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        size = .data[["cluster_size"]] * 1
      ),
      show.legend = FALSE
    ) +
    # color scheme
    ggplot2::scale_fill_manual(
      "Pathway",
      values = col_univ() # nolint
    ) +
    ggplot2::scale_size_area(
      max_size = 16
    ) +
    # text labels
    ggrepel::geom_text_repel(
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        label = .data[["name"]]
      ),
      bg.color = "grey10",
      color = "white",
      bg.r = 0.03,
      alpha = 1,
      nudge_x = 0.25,
      nudge_y = 0.1,
      size = 5
    ) +
    # plot theme
    gen_theme() + # nolint
    net_theme(leg = c(0.8, 0.1)) # nolint
  return(np_anno) # nolint
}
