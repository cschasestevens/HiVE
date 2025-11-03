#' Lipid Network Plot
#'
#' Generates a network plot using known lipid metabolism pathways
#' to visualize class-based enrichment statistics. Uses a manually
#' curated collection of pathways as a reference. Lipid class
#' annotations assigned by this reference enable mapping of any
#' untargeted or targeted lipidomics dataset onto the network.
#'
#' @param level Base level to plot (either 0 or 1). Level 1 contains
#' finer class-based separation while 0 (default) plots higher level
#' lipid classes.
#'
#' @return A ggplot2 object plotting the HiVE base network.
#' @examples
#'
#' # ## Base network
#' # hive_base()
#'
#' @export
hive_base <- function(
  level = 0
) {
  #---- Base network ----
  # format edges
  ne <- nedge[nedge[["Level"]] == level, ] # nolint
  # format nodes
  nn <- nnode[nnode[["Level"]] == level, ] # nolint
  # Create network
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
  ## Plot
  np_base <- ggraph::ggraph(npl) + # Base graph
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
    # node points
    ggraph::geom_node_point() +
    # graph edges and attributes
    ggraph::geom_edge_link(
      ggplot2::aes(label = .data[["enzyme.id"]]), # nolint
      label_dodge = ggplot2::unit(2, "mm"),
      angle_calc = "along",
      alpha = 0.5,
      label_alpha = 0.0,
      color = "grey50"
    ) +
    # color scheme
    ggplot2::scale_fill_manual(
      "Pathway",
      values = col_univ() # nolint
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
  return(np_base) # nolint
}
