#' Lipid Network Plot
#'
#' Generates a network plot using known lipid metabolism pathways
#' to visualize class-based enrichment statistics. Uses a manually
#' curated collection of pathways as a reference. Lipid class
#' annotations assigned by this reference enable mapping of any
#' untargeted or targeted lipidomics dataset onto the network.
#'
#' @param type Network type to plot; current lipid networks available are
#' "base" (Default) and "oxylipins", which is a more detailed subnetwork
#' of HiVE base.
#' @param lab_a Alpha value for enzyme names (default is transparent).
#'
#' @return A ggplot2 object plotting the HiVE base network.
#' @examples
#'
#' # ## Base network
#' # hive_base()
#'
#' @export
hive_base <- function(
  type = "base",
  lab_a = 0.0
) {
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
  # Create network
  np <- igraph::graph_from_data_frame(
    ne,
    vertices = nn
  )
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
  ## Plot
  np_base <- ggraph::ggraph(npl) + # Base graph
    # Subclass shading
    ggplot2::geom_sf(
      data = npw,
      ggplot2::aes(
        fill = .data[["pw"]]
      ),
      color = "grey25",
      alpha = 1
    ) +
    # node points
    ggraph::geom_node_point() +
    # graph edges and attributes
    ggraph::geom_edge_arc(
      ggplot2::aes(label = .data[["col.gene"]]), # nolint
      curvature = 0.1,
      label_dodge = ggplot2::unit(2, "mm"),
      arrow = ggplot2::arrow(length = ggplot2::unit(4, "mm"), type = "closed"),
      start_cap = ggraph::circle(3, "mm"),
      end_cap = ggraph::circle(3, "mm"),
      angle_calc = "along",
      alpha = 0.5,
      label_alpha = 0.5,
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
