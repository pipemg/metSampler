#' library(Hitandrun)
#' @author Felipe de Jesus Muñoz González, Francisco Romero, Mirna Vazquez, Patricia Carvajal
#' @references https://github.com/ComunidadBioInfo/rbioc18/blob/master/materials/collaborative%20projects/projects_descriptions/04_metaboloma.pdf
#' @param two dycon matrix
#' @export
#' @examples   


#Vamos a correr el test para todas las lineas

distribution_results<-function(mat1,mat2){
  x=seq_len(nrow(mat1))###
  resultado=lapply(x, function(x) ks.test(mat1[x,],mat2[x,])$p.value)##distribution_test(mat1, mat2)
  #Si es 0 son diferentes.
  y<-ifelse(as.numeric(resultado) < .01, 1, 0)
  # Accesa a distintas e iguales.
  distintas<-which(y == 0)
  iguales<-which(y == 1)
  return(iguales)
}

