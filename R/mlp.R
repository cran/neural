mlp<-function(inp,weigth,dist,neurons,actfns=c(),layer=NaN){

		sigmoid<-function(x) 1/(1+exp(-x));

		tanhip<-function(x) (1-exp(-2*x))/(1+exp(-2*x));

		gauss<-function(x) exp(-(x^2)/2);

		ident<-function(x) x;

		valuate<-function(fnc){
			value<-list();
			value[1]<-list(inp[watch,]);
			for(i in 2:layer){
				ee<-c();
				for(j in 1:neurons[i]){
					e<-0;
					for(k in 1:neurons[i-1]) e<-e+value[[i-1]][k]*weigth[[i-1]][k,j];
					ee<-c(ee,fnc[[i-1]](e+dist[[i]][j]));
				}
				value[i]<-list(ee);
			}
			value[[layer]];
		}

	if ((length(neurons)!=length(actfns)+1)&(length(actfns)>0)) return("The number of activation function must be the same as the number of layer");
	if (length(actfns)!=0){
		fnd<-FALSE;
		for(i in 1:length(actfns)) 
			if ((actfns[i]!=1)&(actfns[i]!=2)&(actfns[i]!=3)&(actfns[i]!=4)) fnd<-TRUE;
		if (fnd) return("The code of the activation functions must be integer and must be between 1-4");
	}
	if ((is.na(layer))|(layer>=length(neurons))) layer<-length(neurons)
	if (layer<1) layer<-1;
	if (layer==1) return(inp);

	fnc<-c();

	if (length(actfns)>0){
		fnctype<-actfns;
		for(i in 1:(layer-1)){
				if(fnctype[i]==1) {fnc<-c(fnc,sigmoid);}
				if(fnctype[i]==2) {fnc<-c(fnc,tanhip);}
				if(fnctype[i]==3) {fnc<-c(fnc,gauss);}
				if(fnctype[i]==4) {fnc<-c(fnc,ident);}
			}
	}
	else{
		for(i in 1:(layer-1)) fnc<-c(fnc,sigmoid)
		fnctype<-rep(1,times=(layer-1))
	}

	reslt<-c()
	for(watch in 1:nrow(inp))
		reslt<-rbind(reslt,valuate(fnc))
	reslt;
}
