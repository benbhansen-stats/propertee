#' @include Design.R
NULL
# The above ensures that `Design` is defined prior to `DesignStructure`

setClass("DesignStructure",
         contains = "data.frame",
         slots = c(Design = "Design",
                   binary = "logical"))

setValidity("DesignStructure", function(object) {
  validObject(as(object@Design, "Design"))

  # Todo: check validity of `data.frame`. Right dimensions?
  return(TRUE)
})
######### Structure

##' Creates a \code{data.frame} containing the information at the unit of
##' assignment/cluster level.
##' @title Returns \code{Design} Structure Information
##' @param design A \code{Design}.
##' @param binary Binary, default \code{FALSE}. If \code{TRUE} and the design
##'   contains a \code{dichotomy}, replace the treatment column with its binary
##'   representation. Has no effect if \code{design} is not dichotomized.
##' @return A \code{data.frame} containing the structure of the \code{design}.
##' @export
get_structure <- function(design, binary = FALSE) {

  struct <- design@structure
  struct[, design@column_index == "t"] <-
    treatment(design, binary = binary)

  return(new("DesignStructure",
             struct,
             binary = binary,
             Design = design))
}

##' @title Show a \code{DesignStructure}
##' @param object \code{DesignStructure} object
##' @return an invisible copy of \code{object}
##' @export
setMethod("show", "DesignStructure", function(object) {
  struct <- data.frame(object)
  old_row_names <- rownames(struct)

  var_id_row <- object@Design@column_index
  var_id_row <- gsub("^t$", "Treatment", var_id_row)
  uoatype <- switch(object@Design@unit_of_assignment_type,
                    "unit_of_assignment" = "Unit of Assignment",
                    "cluster" = "Cluster",
                    "unitid"  = "Unit ID")
  var_id_row <- gsub("^u$", uoatype, var_id_row)
  var_id_row <- gsub("^b$", "Block", var_id_row)
  var_id_row <- gsub("^f$", "Forcing", var_id_row)

  struct <- rbind(var_id_row, struct)
  rownames(struct) <- c("", old_row_names)


  show(struct)
  invisible(object)
})
