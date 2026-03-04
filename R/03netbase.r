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
#' @param filt_pw Filter by a specific pathway present in either the "base"
#' or "oxylipins" networks.
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
  lab_a = 0.0,
  filt_pw = NULL
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
  npw[["pw"]] <- factor(npw[["pw"]], levels = sort(unique(npw[["pw"]])))
  # Select individual pathway to plot
  if (!is.null(filt_pw)) {
    ## Define colors for each pathway
    col1 <- setNames(
      levels(npw[["pw"]]),
      col_univ()[1:length(levels(npw[["pw"]]))] # nolint
    )
    npw <- npw[npw[["pw"]] == filt_pw, ]
    col1 <- names(col1[grepl(filt_pw, col1)])
  }
  ## Plot
  if (is.null(filt_pw)) {
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
        strength = 0.1,
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
      # color scheme
      ggplot2::scale_fill_manual(
        "Pathway",
        values = col_univ() # nolint
      ) +
      # text labels
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
      # plot theme
      gen_theme() + # nolint
      net_theme(leg = c(0.8, 0.1)) # nolint
  }
  if (!is.null(filt_pw)) {
    np_base <- ggraph::ggraph(npl) + # Base graph
      # Subclass shading
      ggplot2::geom_sf(
        data = npw,
        fill = col1,
        color = "grey25",
        alpha = 1
      ) +
      # node points
      ggraph::geom_node_point(
        ggplot2::aes(
          alpha = ifelse(
            grepl(filt_pw, .data[["fatty.acid"]]),
            1, 0
          )
        )
      ) +
      # graph edges and attributes
      ggraph::geom_edge_arc(
        ggplot2::aes(
          label = .data[["col.gene"]],
          alpha = ifelse(
            grepl(filt_pw, .data[["start.fa"]]),
            1, 0
          )
        ), # nolint
        strength = 0.1,
        label_alpha = lab_a,
        label_dodge = ggplot2::unit(2, "mm"),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(4, "mm"), type = "closed"
        ),
        start_cap = ggraph::circle(3, "mm"),
        end_cap = ggraph::circle(3, "mm"),
        angle_calc = "along",
        color = "grey50"
      ) +
      # color scheme
      ggplot2::scale_fill_manual(
        "Pathway",
        values = col_univ() # nolint
      ) +
      # text labels
      ggplot2::geom_label(
        ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          label = ifelse(
            grepl(filt_pw, .data[["fatty.acid"]]),
            .data[["name"]],
            ""
          ),
          alpha = ifelse(
            grepl(filt_pw, .data[["fatty.acid"]]),
            1, 0
          ),
          linewidth = ifelse(
            grepl(filt_pw, .data[["fatty.acid"]]),
            0.15, NA
          )
        ),
        size = 4,
        color = "white",
        fill = "grey10",
        show.legend = FALSE
      ) +
      # plot theme
      gen_theme() + # nolint
      net_theme(leg = c(0.8, 0.1)) # nolint
  }
  return(np_base) # nolint
}
