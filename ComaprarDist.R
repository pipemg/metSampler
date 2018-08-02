
#Vamos a correr el test para todas las lineas


test <- function(mat1, mat2){
  x=seq_len(nrow(mat1))
  results=lapply(x, function(x) ks.test(mat1[x,],mat2[x,])$p.value)
}

resultado=test(Matriz1, Matriz2)

#Si es 0 son diferentes.


y<-ifelse(as.numeric(resultado) < .01, 1, 0)

# Accesa a distintas e iguales.

distintas<-which(y == 0)
iguales<-which(y == 1)


