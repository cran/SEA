\name{G5Fun}
\alias{G5Fun}
\title{segregation analysis of G5 population}
\description{
Phenotypic observations in G5 population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
G5Fun(df,model,G5text2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{G5text2}{number of plants measured in each family.}
	    }
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G5=data(G5exdata)
G5Fun(G5exdata,"PG-AD",1)
}




