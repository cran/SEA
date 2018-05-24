\name{SEA-package}
\alias{SEA-package}
\alias{SEA}
\docType{package}
\title{
Segregation Analysis
}
\description{
A few major genes and a series of polygenes are responsive for each quantitative trait. Major genes are individually identified while polygenes are collectively detected. This is mixed major genes plus polygenes inheritance analysis or segregation analysis (SEA). In the SEA, phenotypic observations from a single or multiple segregation populations along with their parents are used to fit all the possible models and the best model is viewed as the model of the trait. There are fourteen types of population combinations available. 

}
\details{
\tabular{ll}{
Package: \tab SEA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-5-20\cr
Depends: \tab shiny,MASS,doParallel,foreach\cr
Imports: \tab KScorrect,kolmim,utils,stats,grDevices,graphics,data.table\cr
License: \tab GPL(>=2)\cr
LazyLoad: \tab yes\cr
}
Users can use 'SEA()' start the GUI.
}
\author{
Zhang Ya-Wen, Du Ying-wen,Ren Wen-long and Zhang Yuan-ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
The EIM algorithm in the joint segregation analysis of quantitative traits.
Zhang Yuan-Ming*,Gai Junyi,Yang Yonghua(2003).
}

\keyword{ 
R
SEA
package }
\examples{
\dontrun{ SEA() }
}


