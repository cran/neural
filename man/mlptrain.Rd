\name{mlptrain}
\alias{mlptrain}

\title{MLP neural network}

\description{A simple MLP neural network that is suitable for classification tasks.}

\usage{	mlptrain(inp,neurons,out,alfa=0.2,it=200,online=TRUE,
		permute=TRUE,thresh=0,dthresh=0.1,actfns=c(),visual=TRUE)
}

\arguments{

	\item{inp}{a matrix that contains one input data in each row.}
	\item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{out}{a matrix that contains one output data in each row.}
	\item{alfa}{the learning-rate parameter of the back-propagation algorithm.}
	\item{it}{the maximum number of training iterations.}
	\item{online}{if TRUE the algorithm will operate in sequential mode of back-propagation,if FALSE the algorithm
			will operate in batch mode of back-propagation.}
	\item{permute}{if TRUE the algorithm will use a random
	  permutation of the input data in each epoch.}
	\item{thresh}{the maximal difference between the desired
	  response and the actual response that is regarded as zero.}
	\item{dthresh}{if the difference between the desired response
	  and the actual response is lesser than this value, the
	  corresponding neuron is drawn in red, otherwise it is drawn in green.}
	\item{actfns}{a numeric vector, which contains the numeric code
	  of the activation functions. The length of the vector must be
			the same as the length of the neurons vector, and each element of the vector must be between 1 and 4.
			The possible numeric codes are the following:
			1: Logistic function
			2: Hyperbolic tangent function
			3: Gauss function
			4: Identical function.}
	  	
	\item{visual}{a logical value, that switches on/off the
	  graphical user interface.}
}

\value{
	a list with 4 arguments:

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

}

\details{
	The function creates an MLP neural network on the basis of the
	function parameters. After the creation of the network it is
	trained with the back-propagation algorithm using the inp and
	out parameters. The inp and out parameters has to be the same
	number of rows, otherwise the
	function will stop with an error message. The function has a
	graphical user interface that can be switched on and off using
	the visual argument.
	If the graphical interface is on, the activation functions can
	be set in manually. If the activation functions are not set then
	each of them
	will be automatically the logistic function.
	The result of the function are the parameters of the trained MLP neural network. Use the mlp function for information recall.
}


\seealso{ `mlp' for recall; `rbftrain' and `rbf' for training an RBF network.}

\examples{
	x<-matrix(c(1,1,0,0,1,0,1,0),4,2)
	y<-matrix(c(0,1,1,0),4,1)
	neurons<-4
	\dontrun{
	data<-mlptrain(x,neurons,y,it=4000);
	mlp(x,data$weigth,data$dist,data$neurons,data$actfns)
	}
}

\keyword{neural}
