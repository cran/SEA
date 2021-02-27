\name{G5BCFFun}
\alias{G5BCFFun}
\title{segregation analysis of G5BCF population}
\description{
Phenotypic observations in G5BCF population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
G5BCFFun(df,model,G5BCFtext2)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{G5BCFtext2}{number of plants measured in each family.}
	    }
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G5BCF=data(G5BCFexdata)
G5BCFFun(G5BCFexdata,"1MG-AD",1)
}




