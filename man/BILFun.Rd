\name{BILFun}
\alias{BILFun}
\title{segregation analysis of BIL population}
\description{
Phenotypic observations in BIL population have often been used to identify mixed major-gene plus polygenes inheritance model for quantitative traits in plants.}
\usage{
BILFun(df,model,BILfr)
}
\arguments{
	\item{df}{phenotype matrix.}
	\item{model}{genetic model.}
	\item{BILfr}{BIL type.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
BIL=data(BILexdata)
BILFun(BILexdata,"0MG","BIL1(F1xP1)")
}

