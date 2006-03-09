\name{mlptrain}
\alias{mlptrain}

\title{MLP neural network}

\description{A simple MLP neural network that is suitable for classification tasks.}

\usage{	mlptrain(inp,neurons,out,weigth=c(),dist=c(),alfa=0.2,it=200,online=TRUE,
		permute=TRUE,thresh=0,dthresh=0.1,actfns=c(),diffact=c(),visual=TRUE)
}

\arguments{

	\item{inp}{a matrix that contains one input data in each row.}
	\item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{out}{a matrix that contains one output data in each row.}
	\item{weigth}{the starting weigths of the network.}
	\item{dist}{the starting distortions of the network.}
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
	\item{actfns}{a list, which contains the numeric code
	  of the activation functions, or the activation function itself. The length of the vector must be
			the same as the length of the neurons vector, and each element of the vector must be between 1 and 4 or a function.
			The possible numeric codes are the following:
			1: Logistic function
			2: Hyperbolic tangent function
			3: Gauss function
			4: Identical function.}
	\item{diffact}{a list, which contains the differential of the activation functions. Only need to use if you use your own activation functions.}	
	\item{visual}{a logical value, that switches on/off the
	  graphical user interface.}
}

\value{
	a list with 5 arguments:

	\item{weigth}{the weigths of the network.}
	\item{dist}{the distortions of the network.}
        \item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{actfns}{a list, that contains the activation functions. The length of the list is equal
			 with the number of active layers.}
	\item{diffact}{a list, which contains the differential of the activation functions. The length
			 of the list is equal with the number of active layers.}	
}

\details{
	The function creates an MLP neural network on the basis of the
	function parameters. After the creation of the network it is
	trained with the back-propagation algorithm using the inp and
	out parameters. The inp and out parameters has to be the same
	number of rows, otherwise the function will stop with an error message. 

	If you use the weigth or dist argument, than that variables won't be determined by random. 
        This could be useful if you want to retrain your network. In that case use both of this two arguments. 
        
	From this vesion of the package there is the chance to use your own activation functions, by using the actfns
	argument. If you do this, don't forget to set the differential of the activation functions in the diffact argument at 
	the same order and the same position where you are using the new activation function. (No need of using the diffact argument if you're using
	the preset activation functions.)
 
	The function has a graphical user interface that can be switched on and off using
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
