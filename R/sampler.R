#' Transform the structure to polygon constraints
#' @export
#' 
#' @keywords flux cone
#' 
#' @description 
#' init.sampler gets a matrix and the right-hand side to estabish the constraints of the flux space
#' 
#' @author Muñoz-González, Felipe <fmunoz@lcg.unam.mx>
#' 
#' @Usage init.sampler(mat,rh, transf=NULL)
#' @param mat  Is the matrix that represnets the flux cone
#' @param transf transformation to apply to the matrix 
#' 
#' @return constraint structure  
#' 
#' @examples 
#' A <- rbind(c(1,0), c(0,1), c(-1,0), c(0,-1))
#' B <- c(1, 1, 0, 0)
#' init.sampler(A,B)
#'      

init.sampler<-function(mat, rh, transf=NULL){

  constr <- list(constr=mat, rhs=rh, dir=rep('<=', length(rh)))
  if(is.null(transf))
    transform <- simplex.createTransform(nrow(A))
  constr <- simplex.createConstraints(transform, constr)
  return(constr)
}


#' @example 
#'   A <- rbind(c(1,0), c(0,1), c(-1,0), c(0,-1))
#'   B <- c(1, 1, 0, 0)
#'   cone<-init.sampler(A,B)
#'   run.sampler(cone,1000)  
#'      
run.sampler<-function(cone_constr, num.samples, sampler="hitandrun"){
  if(sampler=="hitandrun"){
    har.init(constr)
    return(hitandrun(constr,num.samples))
  }
}