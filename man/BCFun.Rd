\name{BCFun}
\alias{BCFun}
\title{segregation analysis of BC population}
\description{
Phenotypic observations in BC population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
BCFun(df,model)
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
BC=data(BCexdata)
BCFun(BCexdata,"0MG")
}


