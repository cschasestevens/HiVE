#' Generate Structure from SMILES
#'
#' Creates and/or saves a raster image from a SMILES
#' code using the R CDK framework.
#'
#' @param smi Input SMILES code.
#' @param sav Logical indicating whether the image should be saved
#' to the disk.
#' @param title1 If sav is TRUE, provide a file name for saving
#' the image.
#'
#' @return Raster image generated from SMILEs code.
#' @import rcdk
#' @examples
#'
#'  # view_structure(smi = "C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O")
#'
#' @export
view_structure <- function(
  smi,
  sav = FALSE,
  title1 = NULL
) {
  # Input SMILEs
  mol1 <- rcdk::parse.smiles(smi)
  ## Set plot params
  depictor <- rcdk::get.depictor(width = 800, height = 800, zoom = 1)
  img1 <- rcdk::view.image.2d(mol1[[1]], depictor = depictor)
  ## Plot images
  plot(1:20, 1:20, col = "white", axes = FALSE, ann = FALSE)
  rasterImage(img1, 1, 1, 20, 20)
  # Save image
  if(sav == TRUE) { # nolint
    png(title1, width = 10, height = 10, units = "cm", res = 1000)
    plot(1:20, 1:20, col = "white", axes = FALSE, ann = FALSE)
    par(mar = rep(0, 4))
    rasterImage(img1, 1, 1, 20, 20)
    dev.off()
  }
}
