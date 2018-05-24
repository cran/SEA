\name{F23Fun}
\alias{F23Fun}
\title{segregation analysis of F23 population}
\description{
Phenotypic observations in F23 population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
F23Fun(df,model,m_nf)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{m_nf}{number of plants measured in each family.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
F23=data(F23exdata)
F23Fun(F23exdata,"0MG",1)
}




