#' @export
#' This functions creates a comparation of each flux distribution of two simulations in flux cones
#' and gets the enzymes that creates the same distribution
#' @author Felipe de Jesus Muñoz González
#'  Francisco Romero
#'  Mirna Vazquez
#'  Patricia Carvajal
#' @references https://github.com/ComunidadBioInfo/rbioc18/blob/master/materials/collaborative%20projects/projects_descriptions/04_metaboloma.pdf
#' @param mat1 a m x h matrix that represents the simulations of the flux cone A
#' @param mat2 a m x h matrix that represents the simulations of flux cone B.
#' @examples
#'
#' # distribution_results_iguales(Matriz1, Matriz2)
#'


#Vamos a correr el test para todas las lineas

distribution_results_distintos<-function(mat1,mat2){
  x=seq_len(nrow(mat1))###
  resultado=lapply(x, function(x) ks.test(mat1[x,],mat2[x,])$p.value)##distribution_test(mat1, mat2)
  #Si es 0 son diferentes.
  y<-ifelse(as.numeric(resultado) < .01, 1, 0)
  # Accesa a distintas e iguales.
  distintas<-which(y == 0)
  iguales<-which(y == 1)
  nombres<-row.names(mat1)
  genes_distintos<-nombres[distintas]
  return(genes_distintos)
}



#'
#' @export
#' This functions creates a comparation of each flux distribution of two simulations in flux cones
#' and gets the enzymes that creates different distributions
#' @author Felipe de Jesus Muñoz González
#'  Francisco Romero
#'  Mirna Vazquez
#'  Patricia Carvajal
#' @references https://github.com/ComunidadBioInfo/rbioc18/blob/master/materials/collaborative%20projects/projects_descriptions/04_metaboloma.pdf
#' @param mat1 a m x h matrix that represents the simulations of the flux cone A
#' @param mat2 a m x h matrix that represents the simulations of flux cone B.
#' @examples
#' # distribution_results_distintos(Matriz1, Matriz2)


distribution_results_iguales<-function(mat1,mat2){
  x=seq_len(nrow(mat1))###
  resultado=lapply(x, function(x) ks.test(mat1[x,],mat2[x,])$p.value)##distribution_test(mat1, mat2)
  #Si es 0 son diferentes.
  y<-ifelse(as.numeric(resultado) < .01, 1, 0)
  # Accesa a distintas e iguales.
  distintas<-which(y == 0)
  iguales<-which(y == 1)
  nombres<-row.names(mat1)
  genes_iguales<-nombres[iguales]
  return(genes_iguales)
}
