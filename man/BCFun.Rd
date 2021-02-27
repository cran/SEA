\name{BCFun}
\alias{BCFun}
\title{segregation analysis of BC population}
\description{
Phenotypic observations in BC population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
BCFun(df,model)
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
BC=data(BCexdata)
BCFun(BCexdata,"0MG")
}


