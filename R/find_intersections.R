#' Intersect a Line with Spatial Polygons and Compute Segment Lengths
#'
#' Identifies and retains portions of a line that intersect with spatial polygons,
#' computes the lengths of the intersected segments, and filters out zero-length segments.
#'
#' @param polys An \code{sf} object representing spatial polygons, such as intensity zones or regions.
#' @param line_proj An \code{sf} object representing a projected line geometry, typically the result of transforming a great-circle line between two locations to match the CRS of \code{polys}.
#'
#' @return An \code{sf} object representing the intersected segments, with an added column \code{seg_length} indicating the length of each segment in kilometers.
#'
#' @importFrom sf st_intersection st_length
#' @importFrom units set_units
#' @keywords internal
find_intersections <- function(polys, line_proj){

  ##-- 1. intersect the line with the polygons --------------
  suppressWarnings({
    segments <- st_intersection(polys, line_proj)
  })

  ##-- 2 compute segment lengths ---------------------------
  segments$seg_length <- set_units(st_length(segments), "km")  # or "m"

  ##-- 3 keep only polygons actually hit (length > 0) ------
  segments <- segments[as.numeric(segments$seg_length) > 0, ]

  return(segments)
}
