##' @title Generate Direct Adjusted Weights
##' @param design a Design object created by one of `RCT_Design`, `RD_Design`,
##'   or `Obs_Design`.
##' @param data optionally the data for the analysis to be performed on. May be
##'   excluded if these functions are included as the `weights` argument of a
##'   model.
##' @param unitOfAssignmentIds optional; list connecting names of unit of
##'   assignment variables in `design` to unit of assignment variables in `data`
##' @return a WeightedDesign object
##' @export
##' @rdname WeightCreators
ett <- function(design, data = NULL, unitOfAssignmentIds = NULL) {
  if (is.null(data)) {
    data <- .get_data_from_model(design@call$formula, unitOfAssignmentIds)
  }

  if (!is.null(unitOfAssignmentIds)) {
    design <- .update_unitOfAssignmentIds(design, data, unitOfAssignmentIds)
  }

  .treatment_concordance(design, data)

  #### generate weights

  tx_vec <- design@structure[design@columnIndex == "t"][[1]]
  ind_tx <- levels(tx_vec)[2]

  if(length(levels(tx_vec)) > 2){
    warning("weights for non-binary treatments not yet implemented.")
    weights <- rep(1, nrow(design@structure))
  } else if(!("b" %in% names(table(design@columnIndex)))){
    # If no block is specified, then E_Z is the proportion of clusters who receive
    # the treatment.
    # If a block is specified, then E_Z varies by block and is the proportion
    # of clusters within the block that receive the treatment.
    cluster_df <- data.frame(design@structure[design@columnIndex == "u"],
                             Tx = as.numeric(tx_vec == ind_tx))
    E_Z <- mean(cluster_df$Tx)
    weights <- cluster_df$Tx + (1 - cluster_df$Tx) * E_Z / (1 - E_Z)
  } else{
    # Create a data frame at the block level with the block id, the number of
    # clusters within each block, and the number of clusters receiving the
    # treatment within each block.
    # Then, calculate E_Z.

    block_df <- data.frame(blockid = names(table(design@structure[design@columnIndex == "b"])),
                     block_units = as.numeric(table(design@structure[design@columnIndex == "b"])),
                     tx_units = tapply(as.numeric(tx_vec == ind_tx),
                                       design@structure[design@columnIndex == "b"],
                                       FUN = sum))
    block_df$E_Z <- block_df$tx_units / block_df$block_units

    # Create a cluster-level data frame that merges design structure with block
    # data frame. Add variable Tx that converts treatment to 0/1 for calculation
    # of ETT weights.
    cluster_df <- merge(design@structure, block_df,
                        by.x = colnames(design@structure[which(design@columnIndex == "b")]),
                        by.y = "blockid")
    cluster_df$Tx <- as.numeric(cluster_df[,names(design@structure[design@columnIndex == "t"])] ==
                                  ind_tx)
    weights <- cluster_df$Tx + (1 - cluster_df$Tx) * cluster_df$E_Z / (1 - cluster_df$E_Z)
  }
  .join_design_weights(weights, design, target = "ett", data = data)
}

##' @export
##' @rdname WeightCreators
ate <- function(design, data = NULL, unitOfAssignmentIds = NULL) {
  if (is.null(data)) {
    data <- .get_data_from_model(design@call$formula, unitOfAssignmentIds)
  }

  if (!is.null(unitOfAssignmentIds)) {
    design <- .update_unitOfAssignmentIds(design, data, unitOfAssignmentIds)
  }

  .treatment_concordance(design, data)

  #### generate weights
  tx_vec <- design@structure[design@columnIndex == "t"][[1]]
  ind_tx <- levels(tx_vec)[2]

  if(length(levels(tx_vec)) > 2){
    warning("weights for non-binary treatments not yet implemented.")
    weights <- rep(1, nrow(design@structure))
  } else if(!("b" %in% names(table(design@columnIndex)))){
    # If no block is specified, then E_Z is the proportion of clusters who receive
    # the treatment.
    # If a block is specified, then E_Z varies by block and is the proportion
    # of clusters within the block that receive the treatment.
    cluster_df <- data.frame(design@structure[design@columnIndex == "u"],
                             Tx = as.numeric(tx_vec == ind_tx))
    E_Z <- mean(cluster_df$Tx)
    weights <- cluster_df$Tx / E_Z + (1 - cluster_df$Tx) / (1 - E_Z)
  } else{
    # Create a data frame at the block level with the block id, the number of
    # clusters within each block, and the number of clusters receiving the
    # treatment within each block.
    # Then, calculate E_Z.

    block_df <- data.frame(blockid = names(table(design@structure[design@columnIndex == "b"])),
                           block_units = as.numeric(table(design@structure[design@columnIndex == "b"])),
                           tx_units = tapply(as.numeric(tx_vec == ind_tx),
                                             design@structure[design@columnIndex == "b"],
                                             FUN = sum))
    block_df$E_Z <- block_df$tx_units / block_df$block_units

    # Create a cluster-level data frame that merges design structure with block
    # data frame. Add variable Tx that converts treatment to 0/1 for calculation
    # of ATE weights.
    cluster_df <- merge(design@structure, block_df,
                        by.x = colnames(design@structure[which(design@columnIndex == "b")]),
                        by.y = "blockid")
    cluster_df$Tx <- as.numeric(cluster_df[,names(design@structure[design@columnIndex == "t"])] ==
                                  ind_tx)

    weights <- cluster_df$Tx / cluster_df$E_Z + (1 - cluster_df$Tx) / (1 - cluster_df$E_Z)
  }

  .join_design_weights(weights, design, target = "ate", data = data)
}