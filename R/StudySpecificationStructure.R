#' @include StudySpecification.R
NULL
# The above ensures that `StudySpecification` is defined prior to
# `StudySpecificationStructure`

setClass("StudySpecificationStructure",
         contains = "data.frame",
         slots = c(StudySpecification = "StudySpecification"))

setValidity("StudySpecificationStructure", function(object) {
  validObject(as(object@StudySpecification, "StudySpecification"))

  # Todo: check validity of `data.frame`. Right dimensions?
  return(TRUE)
})
######### Structure

##' @title \code{StudySpecification} Structure Information
##'
##' @description Obtaining a \code{data.frame} which encodes the specification
##'   information.
##' @param specification a \code{StudySpecification} object
##' @return A \code{StudySpecificationStructure} object containing the structure
##'   of the \code{specification} as a \code{data.frame}.
##' @export
##' @rdname StudySpecification_structure
##' @examples
##' data(simdata)
##' spec <- rct_spec(z ~ uoa(uoa1, uoa2) + block(bid), data = simdata)
##' get_structure(spec)
get_structure <- function(specification) {

  struct <- specification@structure
  struct[, specification@column_index == "t"] <-
    treatment(specification)

  return(new("StudySpecificationStructure",
             struct,
             StudySpecification = specification))
}

##' @param object a \code{StudySpecificationStructure} object, typically the
##'   output of \code{get_structure}
##' @export
##' @rdname StudySpecification_structure
setMethod("show", "StudySpecificationStructure", function(object) {
  struct <- data.frame(object)
  old_row_names <- rownames(struct)

  var_id_row <- object@StudySpecification@column_index
  var_id_row <- gsub("^t$", "Treatment", var_id_row)
  uoatype <- switch(object@StudySpecification@unit_of_assignment_type,
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
