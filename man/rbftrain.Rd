\name{rbftrain}
\alias{rbftrain}

\title{RBF neural network}

\description{A simple RBF neural network which suitable for approximation.}

\usage{	rbftrain(inp,neurons,out,weigth=c(),dist=c(),alfa=0.2,it=40,err=0,
        sigma=NaN,online=TRUE,permute=TRUE,visual=TRUE)
}

\arguments{
	\item{inp}{a matrix that contains one input data in each row.}
	\item{neurons}{the number of neurons in the hidden layer.}
	\item{out}{a matrix that contains one output data in each row.}
	\item{weigth}{the starting weigths of the network.}
	\item{dist}{the starting distortions of the network.}
	\item{alfa}{the learning-rate parameter of the back-propagation algorithm.}
	\item{it}{the maximum number of training iterations.}
	\item{err}{the average error at the studying points,if the average error
			anytime lower than this value,the algorithm will stop.}
	\item{sigma}{the width of the Gauss functions.}
	\item{online}{if TRUE the algorithm will operate in sequential mode of back-propagation,if FALSE the algorithm
			will operate in batch mode of back-propagation.}
	\item{permute}{if TRUE the algorithm will use a random
	  permutation of the input data in each epoch.}
	\item{visual}{a logical value, that switches on/off the graphical user interface.}
}

\value{
	list with 4 argument
	\item{weigth}{the weights of the network.}
	\item{dist}{the distortion of the network.}
	\item{neurons}{a numeric vector with length equals to the number of layers in the network, and
			the ith layer will contains neurons[i] neuron.}
	\item{sigma}{the width of the Gauss functions.}
}

\details{
	The function creates an RBF neural network on the basis of the
	function parameters. After the creation of the network the function trains
	it using the back-propagation algorithm using the inp and out parameter. This two parameters row number must be the same, else the
	function will stop with an error message. 

	If you use the weigth or dist argument, than that variables won't be determined by random. 
        This could be useful if you want to retrain your network. In that case use both of this two arguments in the same time. 

	The function works with normalized Gauss-functions, whose width parameter will be the sigma argument. If you want to give the values, this argument 
	should be a matrix, with rows equal the number of neurons in the first layer, and columns equal the number of neurons in the second layer. 
	If the sigma argument is NaN, then the width of each Gauss function will be the half of the distance between the two nearest training samples. If the sigma argument is exactly one number, 
	then all sigma value will be that exact number. 

	The function has a graphical user interface that can be switched
	on and off, with the visual argument. If the graphical user interface is on, then the function could
	show the result of the approximation in a co-ordinate system, if it's a function with one parameter.
 
	The result of the function is the parameters of the trained RBF neural network. Use the rbf function for information recall.
}

\seealso{ `rbf' for recalling; `mlp' and `mlptrain' for classification.}

\examples{

	x<-t(matrix(-5:10*24,1,16));
	y<-t(matrix(sin(pi/180*(-5:10*24)),1,16));
	neurons<-8;
	\dontrun{
	data<-rbftrain(x,neurons,y,sigma=NaN)
	rbf(x,data$weigth,data$dist,data$neurons,data$sigma)
	}
}

\keyword{neural}
