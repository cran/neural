\name{rbf}
\alias{rbf}

\title{RBF neural network}

\description{The recalling method of the RBF network which was trained by the rbftrain function.}

\usage{	rbf(inp,weigth,dist,neurons,sigma)
}

\arguments{

	\item{inp}{a matrix that contains one input data in each row.}
	\item{weigth}{the weigths of the network.}
	\item{dist}{the distortion of the network.}
	\item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{sigma}{the sigma parameters of the network.}
}

\value{ a matrix that contains the response data of the network in each row.}

\details{ the last four argument can be produce by the rbftrain algorithm.}

\seealso{`rbftrain' for training an RBF network, `mlp' and `mlptrain' for classification.}

\keyword{neural}
