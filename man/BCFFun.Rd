\name{BCFFun}
\alias{BCFFun}
\title{segregation analysis of BCF population}
\description{
Phenotypic observations in BCF population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
BCFFun(df,model,BCFtext2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{BCFtext2}{number of plants measured in each family.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
BCF=data(BCFexdata)
BCFFun(BCFexdata,"0MG",1)
}



