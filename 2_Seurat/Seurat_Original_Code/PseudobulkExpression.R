AggregateExpression <- function(
    object,
    assays = NULL,
    features = NULL,
    return.seurat = FALSE,
    group.by = 'ident',
    add.ident = NULL,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    margin = 1,
    verbose = TRUE,
    ...
) {
  return(
    PseudobulkExpression(
      object = object,
      assays = assays,
      features = features,
      return.seurat = return.seurat,
      group.by = group.by,
      add.ident = add.ident,
      layer = 'counts',
      method = 'aggregate',
      normalization.method = normalization.method,
      scale.factor = scale.factor,
      margin = margin,
      verbose = verbose,
      ...
    )
  )
}

#' Averaged feature expression by identity class
#'
#' Returns averaged expression values for each identity class.
#'
#' If layer is set to 'data', this function assumes that the data has been log
#' normalized and therefore feature values are exponentiated prior to averaging
#' so that averaging is done in non-log space. Otherwise, if layer is set to
#' either 'counts' or 'scale.data', no exponentiation is performed prior to  averaging.
#' If \code{return.seurat = TRUE} and layer is not 'scale.data', averaged values
#' are placed in the 'counts' layer of the returned object and 'log1p'
#' is run on the averaged counts and placed in the 'data' layer \code{\link{ScaleData}}
#' is then run on the default assay before returning the object.
#' If \code{return.seurat = TRUE} and layer is 'scale.data', the 'counts' layer contains
#' average counts and 'scale.data' is set to the averaged values of 'scale.data'.
#'
#' @param object Seurat object
#' @param assays Which assays to use. Default is all assays
#' @param features Features to analyze. Default is all features in the assay
#' @param return.seurat Whether to return the data as a Seurat object. Default is FALSE
#' @param group.by Category (or vector of categories) for grouping (e.g, ident, replicate, celltype); 'ident' by default
#' To use multiple categories, specify a vector, such as c('ident', 'replicate', 'celltype')
#' @param add.ident (Deprecated). Place an additional label on each cell prior to pseudobulking
#' @param layer Layer(s) to use; if multiple layers are given, assumed to follow
#' the order of 'assays' (if specified) or object's assays
#' @param slot (Deprecated). Slots(s) to use
#' @param verbose Print messages and show progress bar
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#'
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' If return.seurat is TRUE, returns an object of class \code{\link{Seurat}}.
#' @export
#' @concept utilities
#' @importFrom SeuratObject .FilterObjects
#'
#' @examples
#' data("pbmc_small")
#' head(AverageExpression(object = pbmc_small)$RNA)
#' head(AverageExpression(object = pbmc_small, group.by = c('ident', 'groups'))$RNA)
#'
AverageExpression <- function(
    object,
    assays = NULL,
    features = NULL,
    return.seurat = FALSE,
    group.by = 'ident',
    add.ident = NULL,
    layer = 'data',
    slot = deprecated(),
    verbose = TRUE,
    ...
) {
  return(
    PseudobulkExpression(
      object = object,
      assays = assays,
      features = features,
      return.seurat = return.seurat,
      group.by = group.by,
      add.ident = add.ident,
      layer = layer,
      slot = slot,
      method = 'average',
      verbose = verbose,
      ...
    )
  )
}



## どんな関数になっているか？
PseudobulkExpression
# function (object, ...) 
# {
#   UseMethod(generic = "PseudobulkExpression", object = object)
# }
# <bytecode: 0x14bb801a0>
#   <environment: namespace:Seurat>

methods("PseudobulkExpression")
# [1] PseudobulkExpression.Assay*    PseudobulkExpression.Seurat*   PseudobulkExpression.StdAssay*
#   see '?methods' for accessing help and source code

## ３つのクラスに分かれていることが判明
getS3method("PseudobulkExpression", "Assay")


#' Match the case of character vectors
#'
#' @param search A vector of search terms
#' @param match A vector of characters whose case should be matched
#'
#' @return Values from search present in match with the case of match
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' cd_genes <- c('Cd79b', 'Cd19', 'Cd200')
#' CaseMatch(search = cd_genes, match = rownames(x = pbmc_small))
#'
#'
#'
#' Pseudobulk feature expression by identity class
#'
#' Returns a representative expression value for each identity class
#'
#' @param object Seurat object
#' @param method Whether to 'average' (default) or 'aggregate' expression levels
#' @param assay  The name of the passed assay - used primarily for warning/error messages
#' @param category.matrix A matrix defining groupings for pseudobulk expression 
#' calculations; each column represents an identity class, and each row a sample
#' @param features Features to analyze. Default is all features in the assay
#' @param layer Layer(s) to user; if multiple are given, assumed to follow
#' the order of 'assays' (if specified) or object's assays
#' @param slot (Deprecated) See \code{layer}
#' @param verbose Print messages and show progress bar
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#
#' @return Returns a matrix with genes as rows, identity classes as columns.
#' If return.seurat is TRUE, returns an object of class \code{\link{Seurat}}.
#' @method PseudobulkExpression Assay
#' @rdname PseudobulkExpression
#' @importFrom SeuratObject .IsFutureSeurat
#' @export
#' @concept utilities
#'

# getS3method("PseudobulkExpression", "Assay")
PseudobulkExpression.Assay <- function(
    object,
    assay,
    category.matrix,
    features = NULL,
    layer = 'data',
    slot = deprecated(),
    verbose = TRUE,
    ...
) {
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = '5.0.0',
      what = 'GetAssayData(slot = )',
      with = 'GetAssayData(layer = )'
    )
    layer <- slot
  }
  data.use <- GetAssayData(
    object = object,
    layer = layer
  )
  features.to.avg <- features %||% rownames(x = data.use)
  if (IsMatrixEmpty(x = data.use)) {
    warning(
      "The ", layer, " layer for the ", assay,
      " assay is empty. Skipping assay.", immediate. = TRUE, call. = FALSE)
    return(NULL)
  }
  bad.features <- setdiff(x = features.to.avg, y = rownames(x = data.use))
  if (length(x = bad.features) > 0) {
    warning(
      "The following ", length(x = bad.features),
      " features were not found in the ", assay, " assay: ",
      paste(bad.features, collapse = ", "), call. = FALSE, immediate. = TRUE)
  }
  features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
  if (length(x = features.assay) > 0) {
    data.use <- data.use[features.assay, ]
  } else {
    warning("None of the features specified were found in the ", assay,
            " assay.", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  if (layer == 'data') {
    data.use <- expm1(x = data.use)
    if (any(data.use == Inf)) {
      warning("Exponentiation yielded infinite values. `data` may not be log-normed.")
    }
  }
  data.return <- data.use %*% category.matrix
  return(data.return)
}

#' @method PseudobulkExpression StdAssay
#' @rdname PseudobulkExpression
#' @export
#' @concept utilities
#'

# getS3method("PseudobulkExpression", "StdAssay")
PseudobulkExpression.StdAssay <- function(
    object,
    assay,
    category.matrix,
    features = NULL,
    layer = 'data',
    slot = deprecated(),
    verbose = TRUE,
    ...
) {
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = '5.0.0',
      what = 'GetAssayData(slot = )',
      with = 'GetAssayData(layer = )'
    )
    layer <- slot
  }
  layers.set <- Layers(object = object, search = layer)
  features.to.avg <- features %||% rownames(x = object)
  bad.features <- setdiff(x = features.to.avg, y = rownames(x = object))
  if (length(x = bad.features) > 0) {
    warning(
      "The following ", length(x = bad.features),
      " features were not found in the ", assay, " assay: ",
      paste(bad.features, collapse = ", "), call. = FALSE, immediate. = TRUE)
  }
  features.assay <- Reduce(
    f = intersect,
    x =  c(list(features.to.avg),
           lapply(X = layers.set, FUN = function(l) rownames(object[l]))
    )
  )
  if (length(x = features.assay) == 0) {
    warning("None of the features specified were found in the ", assay,
            " assay.", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  data.return <- as.sparse(
    x = matrix(
      data = 0,
      nrow = length(x = features.assay),
      ncol = ncol(x = category.matrix)
    )
  )
  for (i in seq_along(layers.set)) {
    data.i <- LayerData(object = object,
                        layer = layers.set[i],
                        features = features.assay
    )
    if (layers.set[i] == "data") {
      data.use.i <- expm1(x = data.i)
      if (any(data.use.i == Inf)) {
        warning("Exponentiation yielded infinite values. `data` may not be log-normed.")
      }
    } else {
      data.use.i <- data.i
    }
    category.matrix.i <- category.matrix[colnames(x = data.i),]
    if (inherits(x = data.i, what = 'DelayedArray')) {
      stop("PseudobulkExpression does not support DelayedArray objects")
    } else {
      data.return.i <- as.sparse(x = data.use.i %*% category.matrix.i)
    }
    data.return <- data.return + data.return.i
  }
  return(data.return)
}

#' @param assays Which assays to use. Default is all assays
#' @param return.seurat Whether to return the data as a Seurat object. Default is FALSE
#' @param group.by Categories for grouping (e.g, "ident", "replicate", 
#' "celltype"); "ident" by default
#' @param add.ident (Deprecated) See group.by
#' @param method The method used for calculating pseudobulk expression; one of: 
#' "average" or "aggregate"
#' @param normalization.method Method for normalization, see \code{\link{NormalizeData}}
#' @param scale.factor Scale factor for normalization, see \code{\link{NormalizeData}}
#' @param margin Margin to perform CLR normalization, see \code{\link{NormalizeData}}
#'
#' @method PseudobulkExpression Seurat
#' @rdname PseudobulkExpression
#' @export
#' @concept utilities
#'



# getS3method("PseudobulkExpression", "Seurat")
PseudobulkExpression.Seurat <- function(
    object,
    assays = NULL,
    features = NULL,
    return.seurat = FALSE,
    group.by = 'ident',
    add.ident = NULL,
    layer = 'data',
    slot = deprecated(),
    method = 'average',
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    margin = 1,
    verbose = TRUE,
    ...
) {
  CheckDots(..., fxns = 'CreateSeuratObject')
  if (!is.null(x = add.ident)) {
    .Deprecated(msg = "'add.ident' is a deprecated argument. Please see documentation to see how to pass a vector to the 'group.by' argument to specify multiple grouping variables")
    group.by <- c('ident', add.ident)
  }
  if (!(method %in% c('average', 'aggregate'))) {
    stop("'method' must be either 'average' or 'aggregate'")
  }
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = '5.0.0',
      what = 'AverageExpression(slot = )',
      with = 'AverageExpression(layer = )'
    )
    layer <- slot
  }
  
  if (method == "average") {
    inform(
      message = "As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.",
      .frequency = "once",
      .frequency_id = "AverageExpression"
    )
  }
  
  object.assays <- .FilterObjects(object = object, classes.keep = c('Assay', 'Assay5'))
  assays <- assays %||% object.assays
  
  # `features` is expected to be a 2D array - one vector per assay
  # in the case the `features` a vector, duplicate it for each assay
  if (!inherits(features, what = "list")) {
    features <- rep(list(features), times = length(assays))
  }
  
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(x = assays) == 0) {
      stop("None of the requested assays are present in the object")
    } else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (length(x = layer) == 1) {
    layer <- rep_len(x = layer, length.out = length(x = assays))
  } else if (length(x = layer) != length(x = assays)) {
    stop("Number of layers provided does not match number of assays")
  }
  data <- FetchData(object = object, vars = rev(x = group.by))
  #only keep meta-data columns that are in object
  group.by <- intersect(group.by, colnames(data))
  data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
  if (nrow(x = data) < ncol(x = object)) {
    inform("Removing cells with NA for 1 or more grouping variables")
    object <- subset(x = object, cells = rownames(x = data))
  }
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(x = data),
    FUN = function(i) {
      length(x = levels(x = data[, i]))
    }
  )
  if (any(num.levels == 1)) {
    message(
      paste0(
        "The following grouping variables have 1 value and will be ignored: ",
        paste0(colnames(x = data)[which(num.levels <= 1)], collapse = ", ")
      )
    )
    group.by <- colnames(x = data)[which(num.levels > 1)]
    data <- data[, which(num.levels > 1), drop = F]
  }
  category.matrix <- CreateCategoryMatrix(labels = data, method = method)
  #check if column names are numeric
  col.names <- colnames(category.matrix)
  if (any(!(grepl("^[a-zA-Z]|^\\.[^0-9]", col.names)))) {
    col.names <- ifelse(
      !(grepl("^[a-zA-Z]|^\\.[^0-9]", col.names)),
      paste0("g", col.names),
      col.names
    )
    colnames(category.matrix) <- col.names
    inform(
      message = paste0("First group.by variable `", group.by[1],
                       "` starts with a number, appending `g` to ensure valid variable names"),
      .frequency = "regularly",
      .frequency_id = "PseudobulkExpression"
    )
  }
  
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.return[[assays[i]]] <- PseudobulkExpression(
      object = object[[assays[i]]],
      assay = assays[i],
      category.matrix = category.matrix,
      features = features[[i]],
      layer = layer[i],
      verbose = verbose,
      ...
    )
  }
  if (return.seurat) {
    op <- options(Seurat.object.assay.version = "v5", Seurat.object.assay.calcn = FALSE)
    on.exit(expr = options(op), add = TRUE)
    if (layer[1] == 'scale.data') {
      na.matrix <- as.matrix(x = data.return[[1]])
      na.matrix[1:length(x = na.matrix)] <- NA
      #sum up counts to make seurat object
      summed.counts <- PseudobulkExpression(
        object = object[[assays[1]]],
        assay = assays[1],
        category.matrix = category.matrix,
        features = features[[1]],
        layer = "counts"
      )
      toRet <- CreateSeuratObject(
        counts = summed.counts,
        project = if (method == "average") "Average" else "Aggregate",
        assay = names(x = data.return)[1],
        ...
      )
      LayerData(
        object = toRet,
        layer = "scale.data",
        assay = names(x = data.return)[i]
      ) <- data.return[[1]]
    } else {
      toRet <- CreateSeuratObject(
        counts = data.return[[1]],
        project = if (method == "average") "Average" else "Aggregate",
        assay = names(x = data.return)[1],
        ...
      )
      if (method == "aggregate") {
        LayerData(
          object = toRet,
          layer = "data",
          assay = names(x = data.return)[1]
        ) <- NormalizeData(
          as.matrix(x = data.return[[1]]),
          normalization.method = normalization.method,
          verbose = verbose
        )
      }
      else {
        LayerData(object = toRet,
                  layer = "data",
                  assay = names(x = data.return)[1]
        ) <- log1p(x = as.matrix(x = data.return[[1]]))
      }
    }
    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        if (layer[i] == 'scale.data') {
          summed.counts <- PseudobulkExpression(
            object = object[[assays[i]]],
            assay = assays[i],
            category.matrix = category.matrix,
            features = features[[i]],
            layer = "counts"
          )
          toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = summed.counts)
          LayerData(
            object = toRet,
            layer = "scale.data",
            assay = names(x = data.return)[i]
          ) <- data.return[[i]]
        } else {
          toRet[[names(x = data.return)[i]]] <- CreateAssayObject(
            counts = data.return[[i]],
            check.matrix = FALSE
          )
          if (method == "aggregate") {
            LayerData(
              object = toRet,
              layer = "data",
              assay = names(x = data.return)[i]
            ) <- NormalizeData(
              as.matrix(x = data.return[[i]]),
              normalization.method = normalization.method,
              scale.factor = scale.factor,
              margin = margin,
              verbose = verbose
            )
          }
          else {
            LayerData(
              object = toRet,
              layer = "data",
              assay = names(x = data.return)[i]
            ) <- log1p(x = as.matrix(x = data.return[[i]]))
          }
        }
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
      if (layer[which(DefaultAssay(object = object) %in% names(x = data.return))[1]] != 'scale.data') {
        toRet <- ScaleData(object = toRet, verbose = verbose)
      }
    }
    #add meta-data based on group.by variables
    cells <- Cells(toRet)
    for (i in 1:length(group.by)) {
      if (group.by[i] != "ident") {
        v <- sapply(
          strsplit(cells, "_"),
          function(x) {return(x[i])}
        )
        names(v) <- cells
        toRet <- AddMetaData(toRet,
                             metadata = v,
                             col.name = group.by[i]
        )
      }
    }
    #set idents to pseudobulk variables
    Idents(toRet) <- cells
    
    #make orig.ident variable
    #orig.ident = ident if group.by includes `ident`
    #if not, orig.ident is equal to pseudobulk cell names
    if(any(group.by == "ident")) {
      i = which(group.by == "ident")
      v <- sapply(
        strsplit(cells, "_"),
        function(x) {return(x[i])}
      )
      names(v) <- cells
      toRet <- AddMetaData(toRet,
                           metadata = v,
                           col.name = "orig.ident"
      )
    } else {
      toRet$orig.ident <- cells
    }
    return(toRet)
  } else {
    return(data.return)
  }
}

