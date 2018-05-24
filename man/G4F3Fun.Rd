\name{G4F3Fun}
\alias{G4F3Fun}
\title{segregation analysis of G4F3 population}
\description{
Phenotypic observations in G4F3 population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
G4F3Fun(df,model,G4F3text2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{G4F3text2}{number of plants measured in each family.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G4F3=data(G4F3exdata)
G4F3Fun(G4F3exdata,"PG-AD",1)
}




