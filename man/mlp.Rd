\name{mlp}
\alias{mlp}

\title{MLP neural network}

\description{The recalling method of the MLP network which was trained by the mlptrain function.}

\usage{	mlp(inp,weight,dist,neurons,actfns=c(),layer=NaN, ...) }

\arguments{
	\item{inp}{a matrix that contains one input data in each row.}
	\item{weight}{the weights of the network.}
	\item{dist}{the distortions of the network.}
	\item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{actfns}{a list, which contains the numeric code
	  of the activation functions, or the activation function itself. The length of the vector must be
			the same as the length of the neurons vector, and each element of the vector must be between 1 and 4 or a function.
			The possible numeric codes are the following:
			1: Logistic function
			2: Hyperbolic tangent function
			3: Gauss function
			4: Identical function.}
	\item{layer}{the number of the layer as far as the function is calculating. If the value is NaN, the function will calculate till the output layer.
		     It can be useful if you need the output of one of the hidden layers.}
	\item{\dots}{currently not used.}
}

\value{ a matrix that contains the response data of the network, each row contains one response.}

\details{ the "weight, dist, neurons, actfns" arguments can be determined by the mlptrain algorithm.}

\seealso{`mlptrain' for training an MLP network, `rbf' and `rbftrain' for approximation.}


\keyword{neural}
