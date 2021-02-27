\name{BILFun}
\alias{BILFun}
\title{segregation analysis of BIL population}
\description{
Phenotypic observations in BIL population have often been used to identify mixed major-gene plus polygene inheritance model for quantitative traits in plants.}
\usage{
BILFun(df,model,BILfr)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{BILfr}{BIL type.}
	    }
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
BIL=data(BILexdata)
BILFun(BILexdata,"0MG","BIL1(F1xP1)")
}

