#' Lipid Network Nodes
#'
#' Node metadata for constructing a HiVE lipid network.
#'
#' @format ## `nnode`
#' A data frame with 64 rows and 7 columns:
#' \describe{
#'   \item{Node}{shorthand node name}
#'   \item{Name}{full node name}
#'   \item{synthesis.pathway}{lipid synthesis pathway information}
#'   \item{SetID}{lipid set identity}
#'   \item{Level}{network level}
#'   \item{hive_x}{node x coordinate}
#'   \item{hive_y}{node y coordinate}
#'   ...
#' }
#' @source HiVE v1.10
"nnode"

#' Lipid Network Edges
#'
#' Edge metadata for constructing a HiVE lipid network.
#'
#' @format ## `nedge`
#' A data frame with 68 rows and 8 columns:
#' \describe{
#'   \item{from}{lipid precursor}
#'   \item{to}{lipid product}
#'   \item{synthesis.pathway}{lipid synthesis pathway information}
#'   \item{leading.name}{representative gene/enzyme name}
#'   \item{enzyme.id}{unique gene/enzyme identity code}
#'   \item{hgnc.symbols}{associated genes}
#'   \item{SetID}{lipid set identity}
#'   \item{Level}{network level}
#'   ...
#' }
#' @source HiVE v1.10
"nedge"

#' Example Network Data 1
#'
#' Example data for constructing HiVE network plots
#' based on lipidomics data.
#'
#' @format ## `d`
#' A data frame with 53 rows and 9 columns:
#' \describe{
#'   \item{Cluster.ID}{lipid class ID}
#'   \item{label.saturation}{lipid class stratified by saturation}
#'   \item{p.value}{raw p-value}
#'   \item{FDR}{FDR-adjusted p-value}
#'   \item{n}{cluster size}
#'   \item{Increased}{number of increased lipids within cluster}
#'   \item{Decreased}{number of decreased lipids within cluster}
#'   \item{Ratio}{proportion of increased/decreased lipids within cluster}
#'   \item{label}{lipid class label}
#'   ...
#' }
#' @source HiVE v1.10
"d"

#' Oxylipin Subnetwork Nodes
#'
#' Node metadata for constructing a HiVE oxylipin subnetwork.
#'
#' @format ## `nodeoxy`
#' A data frame with 52 rows and 10 columns:
#' \describe{
#'   \item{Node}{shorthand node name}
#'   \item{Root}{root node}
#'   \item{Name}{full node name}
#'   \item{Synonyms}{node synonyms}
#'   \item{Form}{fatty acid chain length and saturation}
#'   \item{synthesis.pathway}{lipid synthesis pathway information}
#'   \item{enzyme.class}{major enzyme class}
#'   \item{Level}{network level}
#'   \item{hive_x}{node x coordinate}
#'   \item{hive_y}{node y coordinate}
#'   \item{fatty.acid}{starting fatty acid}
#'   ...
#' }
#' @source HiVE v1.10
"nodeoxy"

#' Oxylipin Subnetwork Edges
#'
#' Edge metadata for constructing a HiVE oxylipin subnetwork.
#'
#' @format ## `edgeoxy`
#' A data frame with 56 rows and 8 columns:
#' \describe{
#'   \item{from}{lipid precursor}
#'   \item{to}{lipid product}
#'   \item{synthesis.pathway}{lipid synthesis pathway information}
#'   \item{id.enzyme}{representative enzyme name}
#'   \item{id.gene}{representative gene name}
#'   \item{hgnc.symbols}{associated genes}
#'   \item{SetID}{lipid set identity}
#'   \item{Level}{network level}
#'   \item{major.enzyme}{primary oxylipin enzyme pathway}
#'   ...
#' }
#' @source HiVE v1.10
"edgeoxy"

#' Lipid Gene List
#'
#' Reference lipid gene list.
#'
#' @format ## `netgene`
#' A data frame with 299 rows and 5 columns:
#' \describe{
#'   \item{Gene}{gene name}
#'   \item{SetID}{lipid pathway name}
#'   \item{abbv}{unique pathway code}
#'   \item{RefType}{Reference type}
#'   \item{Source}{DOI/PMCID if applicable}
#'   ...
#' }
#' @source HiVE v1.10
"netgene"