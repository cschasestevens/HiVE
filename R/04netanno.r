#' Lipid Network Annotation Plot
#'
#' Extension of hive_base() for adding annotation count information
#' to a base HiVE network. Node sizes are proportional to the number
#' of individual lipids present in the dataset for each class. See
#' HiVE example data for more details.
#'
#' @param type Network type to plot; current lipid networks available are
#' "base" (Default) and "oxylipins", which is a more detailed subnetwork
#' of HiVE base.
#' @param lab_a Alpha value for enzyme names (default is transparent).
#' @param enr_obj An enrichment object including lipid class labels
#' and the number of lipids within each class. The labels must match
#' nodes from the HiVE base network nodes to properly function. Example
#' data are included for formatting help. Univariate statistical results
#' can also be provided as input for annotating subnetwork plots.
#'
#' @return A ggplot2 object plotting the HiVE annotation network.
#' @examples
#'
#' # ## Annotation count network
#' # hive_anno(enr_obj = d)
#'
#' @export
hive_anno <- function(
  type = "base",
  lab_a = 0.5,
  enr_obj
) {
  if (type == "base") {
    # format edges
    ne <- nedge # nolint
    # format nodes
    nn <- nnode # nolint
    ## gene/enzyme column
    ne[["col.gene"]] <- ne[["id.gene"]]
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
  #---- Annotation network ----
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
  ## Adjust additional node attributes
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
  # Plot Annotations
  ## Plot
  np_anno <- ggraph::ggraph(npl) + # Base graph
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
      strength = 0.1,
      label_dodge = ggplot2::unit(2, "mm"),
      arrow = ggplot2::arrow(length = ggplot2::unit(4, "mm"), type = "closed"),
      start_cap = ggraph::circle(3, "mm"),
      end_cap = ggraph::circle(3, "mm"),
      angle_calc = "along",
      alpha = 0.5,
      label_alpha = lab_a,
      color = "grey50"
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
