context("sanitizeMarkers")

# FIXME Need to check raw `FindMarkers()` and `FindAllMarkers()` return

test_that("seurat_small", {
    expect_message(
        sanitizeMarkers(
            object = seurat_small,
            markers = all_markers_small
        ),
        "Markers are already sanitized"
    )
})
