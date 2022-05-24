# returns the mothers id.
# 
# Args:
#   pedigree: the pedigree object.
get_momid <- function(peddat)
  with(peddat, 
       vapply(mindex, function(x) if(x > 0L) id[x] else NA_integer_, 1L))

# returns the fathers id.
# 
# Args:
#   pedigree: the pedigree object.
get_dadid <- function(peddat)
  with(peddat, 
       vapply(findex, function(x) if(x > 0L) id[x] else NA_integer_, 1L))

# creates an edge list to pass to igraph. An edge is included between children 
# and parents.
# 
# Args:
#   pedigree: the pedigree object.
create_igraph_input <- function(peddat){
  id <- peddat$id
  father <- get_dadid(peddat)
  mother <- get_momid(peddat)
  
  # TODO: this is O(n^2)
  stopifnot(anyDuplicated(id) < 1)
  out <- lapply(id, function(x){
    # find the children
    children_idx <- which(x == father | x == mother)
    children <- if(length(children_idx) > 0)
      id[children_idx] else NULL
    
    # get the correct order (smallest first) and return
    is_larger <- x > children
    
    cbind(
      ifelse(is_larger, children, x        ), 
      ifelse(is_larger, x       , children))
  })
  
  out <- do.call(rbind, out)
  as.data.frame(out[!duplicated(out), ])
}
