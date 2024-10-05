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

#' Regroup idents based on meta.data info
#'
#' For cells in each ident, set a new identity based on the most common value
#' of a specified metadata column.
#'
#' @param object Seurat object
#' @param metadata Name of metadata column
#' @return A Seurat object with the active idents regrouped
#'
#' @export
#' @concept utilities
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- RegroupIdents(pbmc_small, metadata = "groups")
#'