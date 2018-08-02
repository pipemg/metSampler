#Vamos a correr el test para todas las lineas
load("sampler_tables.RData")
row.names(Matriz1)
rownames(Matriz1) <- lapply(as.character(seq_len(200)),paste,"flux",sep="-")
rownames(Matriz2) <- lapply(as.character(seq_len(200)),paste,"flux",sep="-")

test <- function(mat1, mat2){
  x=seq_len(nrow(mat1))
  results=lapply(x, function(x) ks.test(mat1[x,],mat2[x,])$p.value)
}

resultado=test(Matriz1, Matriz2)
head(Matriz1)

#Si es 0 son diferentes

y<-ifelse(as.numeric(resultado) < .01, 1, 0)

# Accesa a distintas e iguales.

distintas <-which(y == 0)
iguales<-which(y == 1)

boxplot(Matriz1[2, ], Matriz2[2, ], Matriz1[1, ], Matriz2[1, ])
boxplot(Matriz1[1, ], Matriz2[1, ])


boxplot(Matriz1[1, ], Matriz2[1, ], names = )

       