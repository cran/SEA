\name{G3DHFun}
\alias{G3DHFun}
\title{segregation analysis of G3DH population}
\description{
Phenotypic observations in G3DH population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
G3DHFun(df,model,G3DHtext2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{G3DHtext2}{number of plants measured in each family.}
	    }
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G3DH=data(G3DHexdata)
G3DHFun(G3DHexdata,"0MG",1)
}




