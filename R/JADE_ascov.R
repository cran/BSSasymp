ASCOV_JADE <- function(sdf,supp=NULL,A=NULL,...)
{
  p<-length(sdf)
  moment4<-NULL   
  moment6<-NULL   
  if(is.null(supp)) supp<-matrix(c(rep(-Inf,p),rep(Inf,p)),ncol=2)
  if(is.null(A)) A<-diag(p)
  
  for(j in 1:p){ 
    moment4[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*x^4}),
                            supp[j,1],supp[j,2],...)$value
    moment6[j]<-integrate(Vectorize(function(x){sdf[[j]](x)*x^6}),
                            supp[j,1],supp[j,2],...)$value   
  }   

  P<-matrix(0,p,p)
  ord<-order(moment4,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]]<-1
  }  
  
  moment4<-moment4[ord]
  moment6<-moment6[ord]

  lambda<-NULL
  for(i in 1:p){
    lambda[i]<-moment4[i]-3
  }  


  ASCOV <- matrix(0,p^2,p^2)
  for(i in 1:p){
   for(j in 1:p){
    if(i!=j){ 
      ASVij<-lambda[i]^2*moment6[i]+lambda[j]^2*moment6[j]+(3*lambda[j]-3*lambda[i]-
             lambda[i]^2)^2-
              
             2*lambda[i]*lambda[j]*moment4[i]*moment4[j]+2*lambda[i]*(3*lambda[j]- 
             3*lambda[i]-lambda[i]^2)*moment4[i]-2*lambda[j]*(3*lambda[j]-3*lambda[i]-
             lambda[i]^2)*moment4[j]
      
      ASVij<-ASVij/((lambda[i]^2+lambda[j]^2)^2)
        
      ASCOVij<-2*lambda[i]*lambda[j]*moment4[i]*moment4[j]-lambda[i]^2*moment6[i]-
               lambda[j]^2*moment6[j]+lambda[i]*(6*lambda[i]-6*lambda[j]+ lambda[i]^2-  
               lambda[j]^2)*moment4[i]+lambda[j]*(6*lambda[j]-6*lambda[i]+lambda[j]^2-
               lambda[i]^2)*moment4[j]+(3*lambda[j]-3*lambda[i]-lambda[i]^2)*
               (3*lambda[i]-3*lambda[j]-lambda[j]^2) 
      
      ASCOVij<-ASCOVij/((lambda[i]^2+lambda[j]^2)^2)
    
      ASCOV<-ASCOV+ASCOVij*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(
             diag(p)[,j],diag(p)[,i]))+ASVij*kronecker(tcrossprod(diag(p)[,j],
             diag(p)[,j]),tcrossprod(diag(p)[,i],diag(p)[,i]))

   }  
     
   if(i==j) ASCOV <- ASCOV+0.25*(moment4[j]-1)*kronecker(tcrossprod(diag(p)[,i],
                diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))    
   } 
  }
    
  W<-crossprod(t(P),solve(A))
  W<-crossprod(diag(sign(rowMeans(W))),W)
  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))
  
  list(W=W,COV_W=COV_W,A=A,COV_A=COV_A)
}

ASCOV_JADE_est <- function(X,mixed=TRUE)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  if(mixed){
    W<-JADE(X)$W
  }else W<-diag(p)
  
  X<-tcrossprod(sweep(X,2,colMeans(X)),W)
 
  moments<-matrix(0,7,p)   
  
  moments[1,]<-1
  moments[3,]<-1

  for(j in 1:p){
    moments[5,j]<-mean(X[,j]^4)
    moments[7,j]<-mean(X[,j]^6)
  } 

  kurt<-rep(0,p)
  for(j in 1:p){
     kurt[j]<-moments[5,j]-3
   }   

  P<-matrix(0,p,p)
  ord<-order(kurt,decreasing=TRUE)
  for(j in 1:p){
    P[j,ord[j]]<-1
  }   
 
  moments<-moments[,ord]
  
  lambda<-rep(0,p)  
  for(i in 1:p){
    lambda[i]<-moments[5,i]-3
  }  

  W<-crossprod(t(P),W)

  ASCOV <- matrix(0,p^2,p^2)
  for(i in 1:p){
   for(j in 1:p){
    if(i!=j){ 
      ASVij<-lambda[i]^2*moments[7,i]+lambda[j]^2*moments[7,j]+(3*lambda[j]-3*lambda[i]-
             lambda[i]^2)^2-
              
             2*lambda[i]*lambda[j]*moments[5,i]*moments[5,j]+2*lambda[i]*(3*lambda[j]- 
             3*lambda[i]-lambda[i]^2)*moments[5,i]-2*lambda[j]*(3*lambda[j]-3*lambda[i]-
             lambda[i]^2)*moments[5,j]
      
      ASVij<-ASVij/((lambda[i]^2+lambda[j]^2)^2)
        
      ASCOVij<-2*lambda[i]*lambda[j]*moments[5,i]*moments[5,j]-lambda[i]^2*moments[7,i]-
               lambda[j]^2*moments[7,j]+lambda[i]*(6*lambda[i]-6*lambda[j]+ lambda[i]^2-  
               lambda[j]^2)*moments[5,i]+lambda[j]*(6*lambda[j]-6*lambda[i]+lambda[j]^2-
               lambda[i]^2)*moments[5,j]+(3*lambda[j]-3*lambda[i]-lambda[i]^2)*
               (3*lambda[i]-3*lambda[j]-lambda[j]^2) 
      
      ASCOVij<-ASCOVij/((lambda[i]^2+lambda[j]^2)^2)
    
      ASCOV<-ASCOV+ASCOVij*kronecker(tcrossprod(diag(p)[,i],diag(p)[,j]),tcrossprod(
             diag(p)[,j],diag(p)[,i]))+ASVij*kronecker(tcrossprod(diag(p)[,j],
             diag(p)[,j]),tcrossprod(diag(p)[,i],diag(p)[,i]))

   }  
     
   if(i==j) ASCOV<-ASCOV+0.25*(moments[5,j]-1)*kronecker(tcrossprod(diag(p)[,i],
                   diag(p)[,i]),tcrossprod(diag(p)[,i],diag(p)[,i]))   
     
   } 
  }
  
  A<-solve(W)
  COV_A<-crossprod(t(tcrossprod(kronecker(diag(p),A),ASCOV)),kronecker(diag(p),t(A)))/n  
  COV_W<-crossprod(t(tcrossprod(kronecker(t(W),diag(p)),ASCOV)),kronecker(W,diag(p)))/n
  
  list(W=W,COV_W=COV_W,A=A,COV_A=COV_A)
}


