#' Lipid Network Flux Plot
#'
#' Extension of hive_base() for plotting ratio analysis results onto
#' a base HiVE network. Solid edges indicate statistically significant
#' ratios, edge colors represent the log2 fold change of the specified
#' comparison, and grey edges indicate that the ratio between two lipid
#' classes was not calculated. See HiVE example data for more details.
#'
#' @param type Network type to plot; current lipid networks available are
#' "base" (Default) and "oxylipins", which is a more detailed subnetwork
#' of HiVE base.
#' @param lab_a Edge label alpha value.
#' @param stat_obj A data frame of statistical results containing
#' fold changes and p-values run by ms_stat_uni using the "HiVE" method.
#' @param col_fc Fold change column name.
#' @param col_p p value column name.
#' @param col_rat Column containing associated lipid class ratios.
#'
#' @return A ggplot2 object plotting a HiVE flux network.
#' @examples
#'
#' # ## Annotation count network
#' # hive_flux(enr_obj = d)
#'
#' @export
hive_flux <- function( # nolint
  type = "base",
  lab_a = 0.5,
  stat_obj,
  col_fc = "Log2FC",
  col_p = "p adj",
  col_rat = "Name"
) {
  if (type == "base") {
    # format edges
    ne <- HiVE::nedge # nolint
    # format nodes
    nn <- HiVE::nnode # nolint
    ## gene/enzyme column
    ne[["col.gene"]] <- ne[["id.gene"]]
    ## pathway column
    nn[["pathway"]] <- nn[["synthesis.pathway"]]
  }
  if (type == "oxylipins") {
    # format edges
    ne <- HiVE::edgeoxy # nolint
    # format nodes
    nn <- HiVE::nodeoxy # nolint
    ## gene/enzyme column
    ne[["col.gene"]] <- ne[["id.gene"]]
    ## pathway column
    nn[["pathway"]] <- nn[["fatty.acid"]]
  }
  # define associated ratio column from stats object
  net1 <- stat_obj
  net1[["associated.ratio"]] <- net1[[col_rat]]
  # format data
  ne <- dplyr::left_join(
    ne,
    net1,
    by = "associated.ratio"
  )
  #---- Lipid ratio network ----
  print("1. Mapping ratio results to HiVE network...")
  ## Adjust edge attributes
  ne[["set_alpha_edge"]] <- ifelse(
    is.na(ne[[col_fc]]),
    1,
    0.25
  )
  ### Scale color from 0-1
  ne[["color_edge"]] <- (
    ne[[col_fc]] - min(ne[[col_fc]], na.rm = TRUE)
  ) / (
    max(ne[[col_fc]], na.rm = TRUE) -
      min(ne[[col_fc]], na.rm = TRUE)
  )
  ne[is.na(ne[["color_edge"]]), ][["color_edge"]] <- 0.5
  ### define color scale
  col_fun <- colorRampPalette(col_grad(scm = 5))
  col_len <- col_fun(length(ne[["color_edge"]]))
  ne[["color_edge_names"]] <- col_len
  ### define linetype for significant ratios
  ne[["ratio_sig"]] <- ifelse(
    is.na(ne[[col_p]]), 1, ne[[col_p]]
  )
  # Plot
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
      ggplot2::aes(,
        alpha = .data[["set_alpha_edge"]]
      ),
      strength = 0.1,
      label_dodge = ggplot2::unit(2, "mm"),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(4, "mm"), type = "closed"
      ),
      start_cap = ggraph::circle(3, "mm"),
      end_cap = ggraph::circle(3, "mm"),
      angle_calc = "along",
      edge_alpha = 0.5,
      label_alpha = lab_a,
      color = "grey50"
    )
  print("4. Add node and edge attributes...")
  # Add nodes based on saturation
  np_enr_nodes <- np_base +
    # node points
    ggraph::geom_node_point(
      ggplot2::aes(alpha = 1),
      color = "grey75",
      show.legend = FALSE
    ) +
    ## Change scale
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    ## Edge attributes showing flux between classes
    ggraph::geom_edge_arc(
      ggplot2::aes(
        alpha = ifelse(
          .data[["color_edge"]] == 0.5,
          0,
          1
        ),
        label = ifelse(
          .data[["color_edge"]] == 0.5,
          "",
          .data[["col.gene"]]
        ),
        color = .data[["color_edge_names"]],
        edge_linetype = ifelse(
          .data[["ratio_sig"]] < 0.05,
          "solid",
          "dashed"
        ),
        edge_width = ifelse(.data[["ratio_sig"]] < 0.05, 1.5, 1)
      ),
      strength = 0.1,
      label_dodge = ggplot2::unit(2, "mm"),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(4, "mm"), type = "closed"
      ),
      start_cap = ggraph::circle(3, "mm"),
      end_cap = ggraph::circle(3, "mm"),
      angle_calc = "along",
      label_alpha = 0.5
    )
  print("5. Add node labels...")
  # Add node labels
  np_enr_labs <- np_enr_nodes +
    ggnewscale::new_scale_color() +
    # ChemRICH node labels
    ## Saturated
    shadowtext::geom_shadowtext(
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        label = .data[["name"]]
      ),
      nudge_x = 0.1,
      nudge_y = 0.1,
      color = "white",
      bg.color = "grey25",
      show.legend = FALSE
    )
  print("6. Add HiVE network theme...")
  # Add theme
  np_enr_thm <- np_enr_labs +
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
  return(np_enr_thm)
}
