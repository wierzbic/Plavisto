
tree_data <- read.table('infiles.csv', sep=",", header=TRUE, stringsAsFactors = FALSE)

maketreelist <- function(df, root = df[1, 1]) {
  #browser()
  r <- list(text = root)
  children = unique(df[df[, 1] == root, 2])
  if(length(children) > 0) {
    r$children <- lapply(children, maketreelist, df = df)
  }
  r
}

tree_data_small <- read.table('infiles.csv', sep=",", header=TRUE, stringsAsFactors = FALSE)
makelist <- function(df){   #Add first level data names; figure out changes
  #r <- list()
  #browser()
  children <- unique(df[,1])
  a <- lapply(children, function(current_level, df){
    #browser()
    b <- list(text = current_level, disabled = TRUE)
    b$children <- lapply(unique(df[df[,1] == current_level, 2]), function(current, df){
      c <- list(text = current)
      c$children <- lapply(unique(df[df[,2] == current, 3]), function(current2) {
        return(list(text = current2))
      })
      return(c)
    }, df = df)
    return(b)
  }, df = df)
  return(a)
}
tmp2 <- (makelist(tree_data_small))


makelist <- function(df){   #Add first level data names; figure out changes
  r <- list(text = "ROOT")
  children <- unique(df[,2])
  r$children <- lapply(children, function(current_level, df){
    r <- list(text = current_level)
    r$children <- lapply(unique(df[df[,2] == current_level, 3]), function(current) {
      return(list(text = current))
    })
    return(r)
  }, df = df)
  return(r)
}

tmp <- (makelist(tree_data_small))

tree_data_small <- tree_data[,1:3]
tmp <- makelist(tree_data_small)










install.packages("data.tree", dependencies = TRUE)
library(data.tree)

tree_data <- read.table('infiles.csv', sep=",", header=TRUE, stringsAsFactors = FALSE)
tree_data_subset <- tree_data[,1:3]
tree_data_subset$pathString <- paste(tree_data_subset$groupname, tree_data_subset$foldername, tree_data_subset$samplename, sep = "/")

tree1 <- as.Node(tree_data_subset)
tree2 <- as.list(as.Node(tree_data_subset), mode = "explicit", uname = FALSE, nameName = "text", keepOnly = "TSS RNAseq Col-0 replicates")


tree_data <- read.table('infiles_tmp.csv', sep=",", header=FALSE, stringsAsFactors = FALSE)

tree_data$pathString <- paste(tree_data$V1, tree_data$V2, sep = "/")
tree2 <- as.list(as.Node(tree_data), mode = "explicit", uname = FALSE, nameName = "text")


nodes <- list(
  list(
    text = "RootA",
    children = list(
      list(
        text = "ChildA1"
      ),
      list(
        text = "ChildA2"
      )
    )
  )
  
)




nodes <- list(
  list(
    text = "RootA",
    children = list(
      list(
        text = "ChildA1"
      ),
      list(
        text = "ChildA2"
      )
    )
  ),
  list(
    text = "RootB",
    children = list(
      list(
        text = "ChildB1"
      ),
      list(
        text = "ChildB2"
      )
    )
  )
)
