\name{BCFFun}
\alias{BCFFun}
\title{segregation analysis of BCF population}
\description{
Phenotypic observations in BCF population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
BCFFun(df,model,BCFtext2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{BCFtext2}{number of plants measured in each family.}
	    }
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
BCF=data(BCFexdata)
BCFFun(BCFexdata,"0MG",1)
}



