# sourceDir sources all of the .R (or .r, .S, .s, .Q, .q) files in the given directory
# sourceDir takes at minimum the path of the directory in which the desired files reside
#' Sources files in the path matching the pattern(s) *.r, *.R, *.S, *.s, *.Q, *.q
#'
#' @param path input string denoting location of files to source
#' @examples
#' sourceDir("/protected/projects/pulmarray/Allegro/COPD_Cancer/functions")
sourceDir <- function(path, trace = TRUE, ...)
{
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}