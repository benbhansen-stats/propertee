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

##' @title \code{Design} Structure Information
##'
##' @description Obtaining a \code{data.frame} which encodes the design
##'   information.
##' @param design a \code{Design} object
##' @param binary default \code{FALSE}. If \code{TRUE} and the design contains a
##'   \code{dichotomy}, replace the treatment column with its binary
##'   representation. Has no effect if \code{design} is not dichotomized.
##' @return A \code{DesignStructure} object containing the structure of the
##'   \code{design} as a \code{data.frame}.
##' @export
##' @rdname Design_structure
##' @examples
##' data(simdata)
##' des <- rct_design(z ~ uoa(cid1, cid2) + block(bid), data = simdata)
##' get_structure(des)
get_structure <- function(design, binary = FALSE) {

  struct <- design@structure
  struct[, design@column_index == "t"] <-
    treatment(design, binary = binary)

  return(new("DesignStructure",
             struct,
             binary = binary,
             Design = design))
}

##' @param object a \code{DesignStructure} object, typically the output of
##'   [get_structure()]
##' @export
##' @rdname Design_structure
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
