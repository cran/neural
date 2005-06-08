mlptrain<-function(inp,neurons,out,alfa=0.2,it=200,online=TRUE,permute=TRUE,thresh=0,dthresh=0.1,actfns=c(),visual=TRUE){
		rect<-function(x,y,value){
			xx<-c(x,x+30,x+30,x)
			yy<-c(y,y,y+30,y+30)
			polygon(xx,yy,col=ifelse(value<=dthresh,"red","green"));
		}

		drawnet<-function(lefut=FALSE,conn=TRUE,fnctype=rep(1,times=length(ls)-1)){
			plot(1:600,1:600,xlab="",ylab="",type="n",axes=FALSE);
			title(main="MLP Network");
			polygon(c(0,65,65,0),c(115,115,75,75),col="lightblue2");
			text(32,99,"Prev",cex=0.9);
			polygon(c(80,145,145,80),c(115,115,75,75),col="lightblue2");
			text(112,99,"Next",cex=0.9);
			polygon(c(540,620,620,540),c(550,550,600,600),col="lightblue2");
			text(580,578,ifelse(lefut==FALSE,"START","EXIT"),cex=0.9);
			text(25,40,"Iteration:",cex=0.9);
			polygon(c(75,160,160,75),c(15,15,55,55),col="turquoise");
			polygon(c(160,240,240,160),c(15,15,55,55),col="snow3");
			text(170,37,ifelse(lefut==FALSE,0,iter),pos=2,cex=0.9);
			text(245,37,it,pos=2,cex=0.9);
			polygon(c(0,80,80,0),c(230,230,190,190),col="ivory");
			text(-10,212,ifelse(conn==FALSE,"CONN"," FNS"),pos=4,cex=0.8);

			if (conn){
				for(i in 1:(length(cordx)-1)){
					for(j in 1:length(cordx[[i]])){
						for(k in 1:length(cordx[[i+1]])){
							xx<-c(cordx[[i]][j]+30,cordx[[i+1]][k]);
							yy<-c(cordy[[i]][j]+15,cordy[[i+1]][k]+15);
							lines(xx,yy,col="blue");
							text((xx[1]+xx[2])/2,(yy[1]+yy[2])/2,round(weigth[[i]][j,k]*100)/100)
						}
					}
				}
			}
			else{
				for(i in 1:(length(ls)-1)){
					cx<-cordx[[i]][length(cordx[[i]])/2+1];
					cy<-cordy[[i]][length(cordy[[i]])/2+1];
					polygon(c(cx+45,cx+115,cx+115,cx+45),c(cy+35,cy+35,cy,cy),col="ivory");
					if(fnctype[i]==1) text(cx+40,cy+17,"SZIGM",pos=4,cex=0.7);
					if(fnctype[i]==2) text(cx+40,cy+17,"TANHIP",pos=4,cex=0.7);
					if(fnctype[i]==3) text(cx+40,cy+17,"EXP",pos=4,cex=0.7);
					if(fnctype[i]==4) text(cx+40,cy+17,"IDENT",pos=4,cex=0.7);
				}
			}

			for(i in 1:length(cordx)){
				for(j in 1:length(cordx[[i]])){
					rect(cordx[[i]][j],cordy[[i]][j],ifelse(i==length(cordx),abs(out[watch,j]-val[[i]][j]),dthresh+1));
					if (i==1)
						text(cordx[[i]][j]-10,cordy[[i]][j]+10,round(val[[i]][j]*100)/100,pos=1)
					else text(cordx[[i]][j]+40,cordy[[i]][j]+10,round(val[[i]][j]*100)/100,pos=1);
				}
			}

		}

		sigmoid<-function(x) 1/(1+exp(-x));

		sigmoiddif<-function(x) (sigmoid(x)*(1-sigmoid(x)));
			##exp(-x)/((1+exp(-x))^2);

		tanhip<-function(x) (1-exp(-2*x))/(1+exp(-2*x));

		tanhipdif<-function(x) 1-tanhip(x)^2;

		gauss<-function(x) exp(-(x^2)/2);

		gaussdif<-function(x) -x*gauss(x);

		ident<-function(x) x;

		identdif<-function(x) 1;
	
		v<-function(l,j) {vv<-0;for(i in 1:ls[l-1]) vv<-vv+val[[l-1]][i]*weigth[[l-1]][i,j];vv+dist[[l]][j]};

		valuate<-function(fnc){
			value<-list();
			value[1]<-list(inp[watch,]);
			for(i in 2:length(ls)){
				ee<-c();
				for(j in 1:ls[i]){
					e<-0;
					for(k in 1:ls[i-1]) e<-e+value[[i-1]][k]*weigth[[i-1]][k,j];
					ee<-c(ee,fnc[[i-1]](e+dist[[i]][j]));
				}
				value[i]<-list(ee);
			}
			value;
		}


		deltaz<-function(fnc,fncdif){
			deltak<-list();
			for(i in length(ls):2){
				dd<-c()
				for(j in 1:ls[i]){
					if (i!=length(ls)) {
						summ<-0;
						for(k in 1:ls[i+1]) summ<-summ+weigth[[i]][j,k]*deltak[[i+1]][k];
						dd<-c(dd,fncdif[[i-1]](v(i,j))*summ);
					}
					else dd<-c(dd,(out[watch,j]-val[[i]][j])*fncdif[[i-1]](v(i,j)));
				}
				deltak[i]<-list(dd);
			}
			deltak;
		}

		permut<-function(v){
			iter<-round(runif(1,min=20,max=40))
			for(i in 1:iter){
				j<-round(runif(1,min=1,max=length(v)))
				k<-round(runif(1,min=1,max=length(v)))
				change<-v[j]
				v[j]<-v[k]
				v[k]<-change
			}
			v;
		}

	if (neurons==0) neurons<-c()

	ls<-c(ncol(inp),neurons,ncol(out))

	if (nrow(inp)!=nrow(out)) return("Different input and output sample number");
	if ((length(ls)!=length(actfns)+1)&(length(actfns)>0)) return("Different activation function and active layer number");
	if (length(actfns)!=0){
		talal<-FALSE;
		for(i in 1:length(actfns)) 
			if ((actfns[i]!=1)&(actfns[i]!=2)&(actfns[i]!=3)&(actfns[i]!=4)) talal<-TRUE;
		if (talal) return("Activation functions is a vector and each element of the vector must be between 1-4.");
	}

	weigth<-list()
	for(i in 1:(length(ls)-1)) {weigth[i]<-list(matrix(c(0),ls[i],ls[i+1]))
					for(j in 1:ls[i])
						for(k in 1:ls[i+1])
							weigth[[i]][j,k]<-runif(1,min=-1,max=1)
				    }

	dist<-list()
	for(i in 2:length(ls)) {k<-c()
				for(j in 1:ls[i]) k<-c(k,runif(1,min=-1,max=1))
				dist[i]<-list(k);
				}

	if (visual){
		cordx<-list();cordy<-list();
		for(i in 1:length(ls)){
			xc<-c();yc<-c();
			for(j in 1:ls[i]){
				xc<-c(xc,300-length(ls)*80+i*160);
				yc<-c(yc,300+ls[i]*30-j*60);
			}
			cordx[i]<-list(xc);cordy[i]<-list(yc);
		}
	}	


	watch<-1;
	fnc<-c();

	if (length(actfns)>0){
		fnctype<-actfns;
		for(i in 1:(length(ls)-1)){
				if(fnctype[i]==1) {fnc<-c(fnc,sigmoid);}
				if(fnctype[i]==2) {fnc<-c(fnc,tanhip);}
				if(fnctype[i]==3) {fnc<-c(fnc,gauss);}
				if(fnctype[i]==4) {fnc<-c(fnc,ident);}
			}
	}
	else{
		for(i in 1:length(ls)-1) fnc<-c(fnc,sigmoid)
		fnctype<-rep(1,times=length(ls)-1)
	}
	val<-valuate(fnc);
	lefut<-FALSE;
	conn<-TRUE;

	if (visual) drawnet(lefut,conn,fnctype);
	
	ext<-FALSE;valt<-FALSE;
	while(!ext){

		if (visual) coor<-locator(1)
		else coor<-list(x=0,y=0);
		if ((!conn)&&(visual)){
			for(i in 1:(length(ls)-1)){
				if((coor$x>cordx[[i]][length(cordx[[i]])/2+1]+45)&(coor$x<cordx[[i]][length(cordx[[i]])/2+1]+115)
				  &(coor$y>cordy[[i]][length(cordy[[i]])/2+1])&(coor$y<cordy[[i]][length(cordy[[i]])/2+1]+35)){
				fnctype[i]<-fnctype[i]%%4+1;
				cx<-cordx[[i]][length(cordx[[i]])/2+1];
				cy<-cordy[[i]][length(cordy[[i]])/2+1];
				polygon(c(cx+45,cx+115,cx+115,cx+45),c(cy+35,cy+35,cy,cy),col="ivory");
				if(fnctype[i]==1) text(cx+40,cy+17,"SZIGM",pos=4,cex=0.7);
				if(fnctype[i]==2) text(cx+40,cy+17,"TANHIP",pos=4,cex=0.7);
				if(fnctype[i]==3) text(cx+40,cy+17,"EXP",pos=4,cex=0.7);
				if(fnctype[i]==4) text(cx+40,cy+17,"IDENT",pos=4,cex=0.7);
				}
			}
		}
		if ((!visual)|(coor$x>540)&(coor$x<600)&(coor$y>550)&(coor$y<600)){
			if (!lefut) {
				fnc<-c()
				fncdif<-c()
				for(i in 1:(length(fnctype))){
					if(fnctype[i]==1) {fnc<-c(fnc,sigmoid);fncdif<-c(fncdif,sigmoiddif);}
					if(fnctype[i]==2) {fnc<-c(fnc,tanhip);fncdif<-c(fncdif,tanhipdif);}
					if(fnctype[i]==3) {fnc<-c(fnc,gauss);fncdif<-c(fncdif,gaussdif);}
					if(fnctype[i]==4) {fnc<-c(fnc,ident);fncdif<-c(fncdif,identdif);}
				}
				
				if (!is.na(it)&(it!=0)){
				iter<-0
				succes<-FALSE
				while((iter<it)&(!succes)){
					iter<-iter+1;
					if (visual){
						polygon(c(75,160,160,75),c(15,15,55,55),col="turquoise");
						text(170,37,iter,pos=2,cex=0.9);
					}
			
					if (permute) perm<-permut(1:nrow(inp)) else perm<-1:nrow(inp);
					w2<-weigth;
					t2<-dist;
					for(ii in 1:nrow(inp)){
						watch<-perm[ii];
						val<-valuate(fnc);

						delta<-deltaz(fnc,fncdif);
						for(k in (length(ls)-1):1){
							for(i in 1:ls[k])
								for(j in 1:ls[k+1]){
									valtoz<-alfa*delta[[k+1]][j]*val[[k]][i];
									w2[[k]][i,j]<-w2[[k]][i,j]+valtoz;
								}			
							for(j in 1:ls[k+1])
								if (online) dist[[k+1]][j]<-dist[[k+1]][j]+alfa*delta[[k+1]][j]
								else t2[[k+1]][j]<-t2[[k+1]][j]+alfa*delta[[k+1]][j];
						}
						if (online) weigth<-w2;
					}
					if (!online){ weigth<-w2;dist<-t2;}
					if (thresh!=0){
						succes<-TRUE;
						watch<-0;
						while((watch<nrow(out))&(succes)){
							watch<-watch+1;
							j<-0;
							val<-valuate(fnc)
							while((j<ncol(out))&(succes)){
								j<-j+1;
								if (abs(out[watch][j]-val[[length(ls)]][j])>thresh) succes<-FALSE
							}
						}
					}
				}
				}
		
				watch<-1;
				lefut<-TRUE;val<-valuate(fnc);
				if (visual) drawnet(lefut,conn,fnctype)
				else ext<-TRUE;
			}
			else ext<-TRUE;
		}

		if (visual){

			if ((coor$x>0)&(coor$x<65)&(coor$y>65)&(coor$y<115)){
				watch<-ifelse(watch>1,watch-1,1);val<-valuate(fnc);drawnet(lefut,conn,fnctype);
			}
			if ((coor$x>80)&(coor$x<145)&(coor$y>65)&(coor$y<115)){
				watch<-ifelse(watch<nrow(inp),watch+1,nrow(inp));val<-valuate(fnc);drawnet(lefut,conn,fnctype);
			}
			if ((coor$x>10)&(coor$x<90)&(coor$y>190)&(coor$y<230)){
				conn<-!conn;drawnet(lefut,conn,fnctype);
			}
			if (valt) {
				cordx[[vi]][vj]<-coor$x;cordy[[vi]][vj]<-coor$y;valt=FALSE;
				drawnet(lefut,conn,fnctype);
			}
			else
			{
			for(i in 1:length(cordx))
				for(j in 1:length(cordx[[i]]))
					if ((coor$x>cordx[[i]][j])&(coor$x<cordx[[i]][j]+30)&(coor$y>cordy[[i]][j])&(coor$y<cordy[[i]][j]+30)){
						valt<-TRUE;
						vi<-i;vj<-j;
					}
			}
		}
	}
	list(weigth=weigth,dist=dist,neurons=ls,actfns=fnctype);
}
