\name{G6FFun}
\alias{G6FFun}
\title{segregation analysis of G6F population}
\description{
Phenotypic observations in G6F population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
G6FFun(df,model,G6Ftext2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{G6Ftext2}{number of plants measured in each family.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G6F=data(G6Fexdata)
G6FFun(G6Fexdata,"PG-AD",1)
}




