\name{G4F2Fun}
\alias{G4F2Fun}
\title{segregation analysis of G4F2 population}
\description{
Phenotypic observations in G4F2 population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
G4F2Fun(df,model)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	    }
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G4F2=data(G4F2exdata)
G4F2Fun(G4F2exdata,"PG-AD")
}




