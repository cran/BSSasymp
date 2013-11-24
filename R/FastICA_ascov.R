ASCOV_FastICAdefl <- function(sdf,gs,dgs,supp=NULL,A=NULL,...)
{
  ng<-length(gs)
  p<-length(sdf)
  if(is.null(supp)) supp<-matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A<-diag(p)
  alphas<-matrix(0,ng,p)
  var_diag<-rep(0,p) 

  for(j in 1:p){
    Ex4<-integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),supp[j,1],supp[j,2])$value
    var_diag[j]<-(Ex4-1)/4

    for(i in 1:ng){
      Eg2<-integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)^2}),
                      supp[j,1],supp[j,2],...)$value
      Eg<-integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)}),
                     supp[j,1],supp[j,2],...)$value
      Egx<-integrate(Vectorize(function(x){sdf[[j]](x)*gs[[i]](x)*x}),
                      supp[j,1],supp[j,2],...)$value
      Edg<-integrate(Vectorize(function(x){sdf[[j]](x)*dgs[[i]](x)}),
                      supp[j,1],supp[j,2],...)$value
     
      alphas[i,j]<-(Eg2-Eg^2-Egx^2)/(Egx-Edg)^2
    }
  }
  
  P<-matrix(0,p,p)
  ba<-rep(0,p)
  for(j in 1:p){
   ba[j]<-min(alphas[,j])
  }
  ord<-order(ba)
  for(j in 1:p){
    P[j,ord[j]]<-1
  }  
  bas<-sort(ba)
  var_diags<-var_diag[ord]
  

  ASV<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASV[i,j]<-bas[i]
      }else if(i==j){
        ASV[i,j]<-var_diags[j]
      }else ASV[i,j]<-bas[j]+1 
    }  
  }

  ASCOV<-diag(as.vector(ASV))

  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASCOV[(i-1)*p+j,(j-1)*p+i]<--bas[i]  
       ASCOV[(j-1)*p+i,(i-1)*p+j]<--bas[i]  
      }else if(i>j){
        ASCOV[(i-1)*p+j,(j-1)*p+i]<--bas[j]
        ASCOV[(j-1)*p+i,(i-1)*p+j]<--bas[j]
      }
    }
  }

  W<-crossprod(t(P),solve(A))
  W<-crossprod(diag(sign(rowMeans(W))),W)
  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W,COV_W=COV_W,A=A,COV_A=COV_A)
}


ASCOV_FastICAdefl_est <- function(X,gs,dgs,mixed=TRUE)
{
  ng<-length(gs)
  n<-dim(X)[1]
  p<-dim(X)[2]
 
  var_diag<-rep(0,p) 
  
  if(mixed){
    fI<-adapt_fICA(X,gs,dgs,rep("",ng))
    W<-fI$W
    alphas<-fI$alphas
  }else{
    alphas<-matrix(0,ng,p)
    for(j in 1:p){
      alphas[,j]<-compute_alphas(X[,j],gs,dgs)
    }
  }  
  
  X<-tcrossprod(sweep(X,2,colMeans(X)),W)   

  for(j in 1:p){
    var_diag[j]<-(mean(X[,j]^4)-1)/4
  }

  P<-matrix(0,p,p)
  ba<-rep(0,p)
  for(j in 1:p){
   ba[j]<-min(alphas[,j])
  }
  ord<-order(ba)
  for(j in 1:p){
    P[j,ord[j]]<-1
  }  
  bas<-sort(ba)
  var_diags<-var_diag[ord]
   
  W<-P%*%W
 
  ASV<-matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASV[i,j]<-bas[i]
      }else if(i==j){
        ASV[i,j]<-var_diags[j]
      }else ASV[i,j]<-bas[j]+1 
    }  
  }

  ASCOV<-diag(as.vector(ASV))

  for(i in 1:p){
    for(j in 1:p){
      if(i<j){ 
        ASCOV[(i-1)*p+j,(j-1)*p+i]<--bas[i]  
        ASCOV[(j-1)*p+i,(i-1)*p+j]<--bas[i]  
      }else if(i>j){
        ASCOV[(i-1)*p+j,(j-1)*p+i]<--bas[j]
        ASCOV[(j-1)*p+i,(i-1)*p+j]<--bas[j]
      }
    }
  }

  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W,COV_W=COV_W,A=A,COV_A=COV_A)
}



compute_alphas<-function(x,gs,dgs)
{
  alphas <- NULL
  for(i in 1:length(gs)){
   Eg <- mean(gs[[i]](x))
   Eg2 <- mean(gs[[i]](x)^2)
   Egx <- mean(gs[[i]](x)*x)
   Edg <- mean(dgs[[i]](x))
   if((Egx-Edg)==0){
    alphas[i] <- Inf
   }else alphas[i] <- (Eg2-Eg^2-Egx^2)/(Egx-Edg)^2
    
   if(alphas[i]<0) alphas[i] <- Inf
  }
  alphas
}

