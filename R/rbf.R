rbf<-function(inp,weigth,dist,neurons,sigma){


		gauss<-function(x,sigma) exp(-(x^2)/(2*sigma^2));

		ident<-function(x) x;


		valuate3<-function(pont){
			value<-list();
			value[1]<-list(pont);
			i<-2;
			total<-0;
			ee<-c();


			for(j in 1:neurons[i]){
				e<-0;
				for(k in 1:neurons[i-1]) 
					e<-e+abs(value[[i-1]][k]-weigth[[i-1]][k,j]);
				total<-total+gauss(e,ifelse((is.nan(sigma))&(neurons[1]==1),dist[[i]][j],sigma));
			}


			for(j in 1:neurons[i]){
				e<-0;
				for(k in 1:neurons[i-1]) 
					e<-e+abs(value[[i-1]][k]-weigth[[i-1]][k,j]);
				ee<-c(ee,gauss(e,ifelse((is.nan(sigma))&(neurons[1]==1),dist[[i]][j],sigma))/total);
			}
			value[i]<-list(ee);

			i<-3
			e<-0;
			for(k in 1:neurons[i-1]) 
				e<-e+value[[i-1]][k]*weigth[[i-1]][k,1];
			ident(e+dist[[i]][1]);
		}



	value<-c();
	for(i in 1:nrow(inp)) value<-rbind(value,valuate3(inp[i,]));
	value;
	
}
