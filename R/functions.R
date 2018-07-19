

#' Title
#'
#' @param fileName 
#'
#' @return
#' @export
#'
#' @examples
readEcoNetwork<- function(fileName)
{
 web <- read.csv(fileName)
  if ((ncol(web)-1)==nrow(web))  
    {
    tmp0 <- ncol(web)
    tmp1 <- web[,2:tmp0]
    tmp2 <- as.matrix(tmp1)
    g <- graph_from_adjacency_matrix(tmp2)
    }  else 
        {
        web <- read.csv(fileName,header = T,check.names = F)
        web <- web[,c(2,1)]
        g <- graph_from_data_frame(web)  
        }
}  




#' Title
#'
#' @param g 
#'
#' @return
#' @export
#'
#' @examples
plotEcoNetworkTrophLevel <- function(g)
{
  
  require(NetIndices)
  
  matadj <- get.adjacency(g,sparse = F)
  
  #matadj
  
  tl <- TrophInd(matadj)#trophic level and omnivory
  
  
# Plot with trophic levels

lMat <-matrix(
  nrow=vcount(g),  # Rows equal to the number of vertices (species)
  ncol=2)
  
  lMat[,2]<-jitter(tl$TL,0.1)   # y-axis value based on trophic level
  lMat[,1]<-runif(vcount(g))   # randomly assign along x-axis
  
#  plot(g, edge.width=.3,edge.arrow.size=.4,
 #      vertex.label=NA,
  #     edge.color="grey50",
    #   edge.curved=0.3, layout=lMat)
  
  require(RColorBrewer)
  
  colTL <- as.numeric(cut(tl$TL,11))#divide trophic levels in 11 segments
  colnet <-  brewer.pal(11,"RdYlGn")#assign colors to trophic levels
  V(g)$color <- colnet[12-colTL]#add colors to the igraph object
  
  plot(g, edge.width=.3,edge.arrow.size=.4,
       vertex.label=NA,
       edge.color="grey50",
       edge.curved=0.3, layout=lMat)  
  
}




#' Title
#'
#' @param g 
#'
#' @return
#' @export
#'
#' @examples
topologicalIndicesEcoNetwork <- function(g)
{
  
  #calculate other network indices
  
  deg <- degree(g, mode="out") # calculate the in-degree: the number of predators  
  
  V(g)$outdegree <-  deg
  
  nTop <- length(V(g)[outdegree==0]) # Top predators do not have predators
  
  deg <- degree(g, mode="in") # calculate the in-degree: the number of preys
  
  V(g)$indegree <-  deg
  
  nBasal <- length(V(g)[indegree==0]) # Basal species do not have preys 
  
  vcount(g)-nTop-nBasal # Intermediate species
  
  size <- vcount(g)
  
  links <- ecount(g)
  
  linkDen <- links/size          # Linkage density
  
  conn <- links/size^2           # Connectance
  
  pathLength <- average.path.length(g)   # Average path length 
  
  clusCoef <- transitivity(g, type = "global")   # Clustering coefficient
  
  #Devolver un data.frame
  
  
  return (data.frame(,size,connectance=conn,linkDen))
}  


