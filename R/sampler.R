#' Iniciate a sampler "metSampler"
#' @export
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

  require(hitandrun)

  constr <- list(constr=mat, rhs=rh, dir=rep('<=', length(rh)))
  if(is.null(transf))
    transform <- simplex.createTransform(nrow(A))
  constr <- simplex.createConstraints(transform, constr)
  return(constr)
}

#' run.sampler "metSampler"
#'
#' @export
#'
#' @description
#' Method that generates n random points in a multidimentional space whose stable state converges on the uniform
#' distribution over a convex polytope defined by a set of linear inequality constraints#'
#' @author Muñoz-González, Felipe <fmunoz@lcg.unam.mx>
#'
#' @Usage run.sampler(cone=NULL,num.samples=NULL,sampler="hitandrun")
#'
#' @param cone_constr  Linear constraints that define the sampling space
#' @param num.samples The desired number of samples to return.
#' @param sampler Selection of the sampler method to use
#'
#' @return constraint structure list of lists
#'
#' @examples
#' A <- rbind(c(1,0), c(0,1), c(-1,0), c(0,-1))
#' B <- c(1, 1, 0, 0)
#' init.sampler(A,B)
#' run.sampler(cone,1000)
#'
#'
run.sampler<-function(cone_constr, num.samples, sampler="hitandrun"){
  require(hitandrun)
  if(sampler=="hitandrun"){
    har.init(constr)
    return(hitandrun(constr,num.samples))
  }
}
