\name{metSampler-package}
\alias{metSampler-package}
\alias{metSampler}
\docType{package}
\title{
  metSampler
}
\description{
  Receives two k-cones from the package Dycon that represent two different null spaces
       from a metabolic constraint based model and generate box plot graphics of the results.
}
\details{
  The package cosists of three modules. The first module, samples two different null spaces from a metabolic contraint based model. This is accomplished through the Hit and Run algorithm (Hit and Run Cran package). The second module tests for diferences in distribution within the two different null spaces for each enzime. This test is performed through the Kolmogorov Smirnov test, the threshold on the p-value for accepting or rejecting the H0 is set to p=0.01. The last modulel generates a box plot for each of the null space enzimes.
}
\author{
Felipe de Jesus Muñoz González, Francisco Romero, Mirna Vazquez, Patricia Carvajal.

Maintainer: Felipe de Jesus Muñoz González <fmunoz@lcg.unam.mx>
}
\references{
  https://github.com/ComunidadBioInfo/rbioc18/blob/master/materials/collaborative%20projects/projects_descriptions/04_metaboloma.pdf.
}
\keywords{ Dycon, Random Sampling, Flux Space of Metabolic Models }
\seealso{
  https://rdrr.io/github/cdiener/dycone/
}
\examples{
distribution_results_distintos(Matriz1, Matriz2)
distribution_results_iguales(Matriz1, Matriz2)

  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
  }
}
