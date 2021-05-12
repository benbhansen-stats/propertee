setClass("Design",
         slots = c(structure = "data.frame",
                   columnIndex = "factor"))

setValidity("Design", function(object) {
  if (any(dim(object@structure) == 0)) {
    "@structure must have positive dimensions"
  } else if (ncol(object@structure) != length(object@columnIndex)) {
    "@columnIndex does not agree with number of columns in @structure"
  } else {
    TRUE
  }
})

RCT_Design <- function(form, data, subset = NULL) {
  if (!is.null(subset)) {
    d <- subset(data, subset = subset)
  }
  m <- as.data.frame(as.matrix(model.frame(form, data)))
  # TODO: Ensure only one `cluster()` and only one `block()` in formula
  tt <- terms(form, c("cluster", "block"))

  index <- factor(rep("t", ncol(m)),
                  levels = c("t", "c", "b", "f"))

  # Handle clusters
  clusters <- grepl("^cluster", colnames(m))
  index[which(clusters)] <- "c"
  cvars <- colnames(m)[clusters][1]
  cvars <- sub("^cluster\\(", "", cvars)
  cvars <- sub("\\)[\\.0-9]*$", "", cvars)
  cvars <- gsub(" ", "", cvars)
  colnames(m)[clusters] <- strsplit(cvars, ",")[[1]]

  # Handle blocks
  blocks <- grepl("^block", colnames(m))
  index[which(blocks)] <- "b"
  bvars <- colnames(m)[blocks][1]
  bvars <- sub("^block\\(", "", bvars)
  bvars <- sub("\\)[\\.0-9]*$", "", bvars)
  bvars <- gsub(" ", "", bvars)
  colnames(m)[blocks] <- strsplit(bvars, ",")[[1]]

  # Handle forcing
  forcings <- grepl("^forcing", colnames(m))
  index[which(forcings)] <- "f"
  fvars <- colnames(m)[forcings][1]
  fvars <- sub("^forcing\\(", "", fvars)
  fvars <- sub("\\)[\\.0-9]*$", "", fvars)
  fvars <- gsub(" ", "", fvars)
  colnames(m)[forcings] <- strsplit(fvars, ",")[[1]]

  new("Design",
      structure = m,
      columnIndex = index)
}

forcing <- block <- cluster <- function (...)
{
  #browser()
  allf <- list(...)
  do.call(cbind, allf)
}
