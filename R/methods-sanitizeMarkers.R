#' Sanitize Markers Output
#'
#' @note [Seurat::FindAllMarkers()] maps the counts matrix rownames correctly in
#'   the `gene` column, whereas [Seurat::FindMarkers()] maps them correctly in
#'   the rownames of the returned marker `data.frame`.
#'
#' @name sanitizeMarkers
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param markers Original [Seurat::FindAllMarkers()] `data.frame`.
#'
#' @return `data.frame`, arranged by adjusted P value.
#'
#' @examples
#' # seurat ====
#' object <- seurat_small
#'
#' # `FindAllMarkers` return
#' invisible(capture.output(
#'     all_markers <- Seurat::FindAllMarkers(object)
#' ))
#' all_sanitized <- sanitizeMarkers(
#'     object = object,
#'     markers = all_markers
#' )
#' glimpse(all_sanitized)
#'
#' # `FindMarkers()` return
#' invisible(capture.output(
#'     ident_3_markers <- Seurat::FindMarkers(
#'         object = object,
#'         ident.1 = "3",
#'         ident.2 = NULL
#'     )
#' ))
#' ident_3_sanitized <- sanitizeMarkers(
#'     object = object,
#'     markers = ident_3_markers
#' )
#' glimpse(ident_3_sanitized)
NULL



# Methods ======================================================================
#' @rdname sanitizeMarkers
#' @export
setMethod(
    "sanitizeMarkers",
    signature("seurat"),
    function(object, markers) {
        validObject(object)
        assert_is_data.frame(markers)

        # Early return on sanitized data
        if (.isSanitizedMarkers(markers, package = "Seurat")) {
            message("Markers are already sanitized")
            return(markers)
        }

        assertHasRownames(markers)
        data <- markers

        # Map the matrix rownames to `rownames` column in tibble
        if ("cluster" %in% colnames(data)) {
            # `FindAllMarkers()` return
            all <- TRUE
            assert_is_subset("gene", colnames(data))
            data <- remove_rownames(data)
            data[["rowname"]] <- data[["gene"]]
            data[["gene"]] <- NULL
        } else {
            # `FindMarkers()` return
            all <- FALSE
            data <- rownames_to_column(data)
        }

        # Now ready to coerce
        data <- data %>%
            as_tibble() %>%
            # Sanitize column names into lowerCamelCase
            camel()

        # Update legacy columns
        if ("avgDiff" %in% colnames(data)) {
            message(paste(
                "Renaming legacy `avgDiff` column to `avgLogFC`",
                "(changed in Seurat v2.1)"
            ))
            data[["avgLogFC"]] <- data[["avgDiff"]]
            data[["avgDiff"]] <- NULL
        }

        # Rename P value columns to match DESeq2 conventions
        if ("pVal" %in% colnames(data)) {
            data[["pvalue"]] <- data[["pVal"]]
            data[["pVal"]] <- NULL
        }
        if ("pValAdj" %in% colnames(data)) {
            data[["padj"]] <- data[["pValAdj"]]
            data[["pValAdj"]] <- NULL
        }

        # Add Ensembl gene IDs
        gene2symbol <- suppressMessages(gene2symbol(object))

        # Only attempt to join rowRanges if we have gene-to-symbol mappings
        if (!is.null(gene2symbol)) {
            gene2symbol[["rowname"]] <- make.unique(gene2symbol[["geneName"]])
            gene2symbol[["geneName"]] <- NULL
            data <- left_join(
                x = data,
                y = gene2symbol,
                by = "rowname"
            )

            # Check that all rows match a geneID
            stopifnot(!any(is.na(data[["geneID"]])))

            # Add rowData
            rowData <- rowData(object)
            rowData <- camel(rowData)
            # Ensure any nested list columns are dropped
            cols <- vapply(
                X = rowData,
                FUN = function(x) {
                    !is.list(x)
                },
                FUN.VALUE = logical(1L)
            )
            rowData <- rowData[, cols]
            rowData <- as.data.frame(rowData)
            data <- left_join(data, rowData, by = "geneID")

            # Ensure that required columns are present
            requiredCols <- c(
                "rowname",
                "geneID",
                "geneName",
                "pct1",
                "pct2",
                "avgLogFC",     # Seurat v2.1
                "padj",
                "pvalue"        # Renamed from `p_val`
            )
            assert_is_subset(requiredCols, colnames(data))
        }

        if (isTRUE(all)) {
            # `cluster` is only present in `FindAllMarkers() return`
            data <- data %>%
                select(!!sym("cluster"), everything()) %>%
                group_by(!!sym("cluster")) %>%
                arrange(!!sym("padj"), .by_group = TRUE)
        } else {
            data <- data %>%
                arrange(!!sym("padj")) %>%
                as.data.frame() %>%
                column_to_rownames()
        }

        data
    }
)
