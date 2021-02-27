\name{G6Fun}
\alias{G6Fun}
\title{segregation analysis of G6 population}
\description{
Phenotypic observations in G6 population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
G6Fun(df,model)
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
G6=data(G6exdata)
G6Fun(G6exdata,"PG-AD")
}




