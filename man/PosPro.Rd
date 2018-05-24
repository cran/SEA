\name{PosPro}
\alias{PosPro}
\title{Posterior Probability}
\description{
calculate posterior probability of the optimal model}
\usage{
PosPro(Population,result,data)}
\arguments{
	\item{Population}{which Population to analysis.}
	\item{result}{result of calculation used corresponding population    function.}
	\item{data}{phenotype matrix.}
	    }
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
F23=data(F23exdata)
result<-F23Fun(F23exdata,"1MG-AD",1)
PosPro("F2:3",result,F23exdata)
}




