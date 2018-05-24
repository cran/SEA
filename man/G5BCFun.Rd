\name{G5BCFun}
\alias{G5BCFun}
\title{segregation analysis of G5BC population}
\description{
Phenotypic observations in G5BC population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
G5BCFun(df,model)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G5BC=data(G5BCexdata)
G5BCFun(G5BCexdata,"1MG-AD")
}




