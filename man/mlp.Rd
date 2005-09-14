\name{mlp}
\alias{mlp}

\title{MLP neural network}

\description{The recalling method of the MLP network which was trained by the mlptrain function.}

\usage{	mlp(inp,weigth,dist,neurons,actfns=c(),layer=NaN)}

\arguments{
	\item{inp}{a matrix that contains one input data in each row.}
	\item{weigth}{the weigths of the network.}
	\item{dist}{the distortions of the network.}
	\item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{actfns}{a numeric vector, that contains the numeric codes
	  of the activation functions. The length of the vector must be
			the same as the length of the neurons vector, and each element of the vector must be between 1-4.
			The possible numeric codes are the following:
			1: Logistic function
			2: Hyperbolic tangent function
			3: Gauss function
			4: Identical function.}
	\item{layer}{the number of the layer as far as the function is calculating. If the value is NaN, the function will calculate till the output layer.
		     It can be useful if you need the output of one of the hidden layers.}
}

\value{ a matrix that contains the response data of the network, each row contains one response.}

\details{ the last four argument is determined by the mlptrain algorithm.}

\seealso{`mlptrain' for training an MLP network, `rbf' and `rbftrain' for approximation.}


\keyword{neural}
