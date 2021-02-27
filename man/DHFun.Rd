\name{DHFun}
\alias{DHFun}
\title{segregation analysis of DH population}
\description{
Phenotypic observations in DH population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
DHFun(df,model)
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
DH=data(DHexdata)
DHFun(DHexdata,"0MG")
}




