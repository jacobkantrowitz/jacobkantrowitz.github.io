# sourceDir sources all of the .R (or .r, .S, .s, .Q, .q) files in the given directory
# sourceDir takes at minimum the path of the directory in which the desired files reside
sourceDir <- function(path, trace = TRUE, ...)
{
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}