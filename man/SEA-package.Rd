\name{SEA-package}
\alias{SEA-package}
\alias{SEA}
\docType{package}
\title{
Segregation Analysis
}
\description{
A few major genes and a series of polygene are responsive for each quantitative trait. Major genes are individually identified while polygene is collectively detected. This is mixed major genes plus polygene inheritance analysis or segregation analysis (SEA). In the SEA, phenotypes from a single or multiple bi-parental segregation populations along with their parents are used to fit all the possible models and the best model for population phenotypic distributions is viewed as the model of the trait. There are fourteen types of population combinations available. Zhang YM, Gai JY, Yang YH (2003) <doi:10.1017/S0016672303006141>.

}
\details{
\tabular{ll}{
Package: \tab SEA\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2021-02-24\cr
Depends: \tab shiny,MASS,doParallel,foreach\cr
Imports: \tab KScorrect,kolmim,utils,stats,grDevices,graphics,data.table\cr
License: \tab GPL(>=2)\cr
LazyLoad: \tab yes\cr
}
Users can use 'SEA()' start the GUI.
}
\author{
Wang Jing-Tian, Zhang Ya-Wen, and Zhang Yuan-Ming \cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
The EIM algorithm in the joint segregation analysis of quantitative traits.
Zhang Yuan-Ming*,Gai Junyi,Yang Yonghua(2003).
}
\examples{
\dontrun{ SEA() }
}


