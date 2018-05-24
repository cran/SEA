\name{F2Fun}
\alias{F2Fun}
\title{segregation analysis of F2 population}
\description{
Phenotypic observations in F2 population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
F2Fun(df,model)
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
F2=data(F2exdata)
F2Fun(F2exdata,"0MG")
}




