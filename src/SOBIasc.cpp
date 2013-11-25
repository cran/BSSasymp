#include <cmath>
//#include <R.h>
//#include <Rdefines.h>
//#include <Rinternals.h>


 extern "C" {

 using namespace std;

 double **prepmat(double *X, int n, int k)
  //be careful to free this memory!!!
  {
    int i;
    int j;
    double **Y = new double* [n];
    for (i=0; i<n; i++) Y[i]=new double [k];
    for (i=0; i<n; i++) 
      for (j=0; j<k; j++)
    Y[i][j]=X[j*n+i];
    return Y;
  }

  int abs(int k)
  {
   int res;
   if(k>0){ 
    res = k;  
   }else res = -k;
   
   return res; 
  }

  double D_lm(double *F, int p, int q, int i, int j, int l, int m, double *B)
  {
    int k;
    int lm;
    if(l>m){
     lm = l;
    }else{ 
     lm = m;
    }
    double D = 0;
    if(i!=j){
     for(k=(-q+lm); k<(q-lm); k++){ 
      D += 0.5*(F[i+i*p+abs(k+l)*p*p]*F[j+j*p+abs(k+m)*p*p]+
           F[i+i*p+abs(k+l)*p*p]*F[j+j*p+abs(k-m)*p*p]);   
     }
     D += 0.25*(B[i+j*p]-1)*(F[i+j*p+l*p*p]+F[j+i*p+l*p*p])*
            (F[i+j*p+m*p*p]+F[j+i*p+m*p*p]); 
    }else{
     for(k=(-q+lm); k<(q-lm); k++){ 
      D += F[i+i*p+abs(k+l)*p*p]*F[i+i*p+abs(k+m)*p*p]+
           F[i+i*p+abs(k+l)*p*p]*F[i+i*p+abs(k-m)*p*p];   
     }
     D += (B[i+i*p]-3)*F[i+i*p+l*p*p]*F[i+i*p+m*p*p]; 
    }
 
    return D;
  }

  
  void ascov(double *F, double *Lambda, double *taus, int *nk, double *B, double *result)
  {
    int i=nk[0];
    int j=nk[1];
    int p=nk[2];
    int q=nk[3];
    int K=nk[4];
    int k;
    int l;
    int m;
 
    double ASV = 0;
    double ASCOV = 0;    
    double sl = 0;
    double sl2 = 0;
  
    for(k=0; k<K; k++){
      ASV += (Lambda[i+i*p+k*p*p]-Lambda[j+j*p+k*p*p])*(Lambda[i+i*p+k*p*p]-
             Lambda[j+j*p+k*p*p])*D_lm(F,p,q,i,j,taus[k],taus[k],B);
    }
 
    if(K>1){
    for(l=0; l<(K-1); l++){
      for(m=(l+1); m<K; m++){
          ASV += 2*(Lambda[i+i*p+l*p*p]-Lambda[j+j*p+l*p*p])*(Lambda[i+i*p+m*p*p]-
                      Lambda[j+j*p+m*p*p])*D_lm(F,p,q,i,j,taus[l],taus[m],B);
        }
       }
      }      
     
     for(k=0; k<K; k++){
      sl += Lambda[i+i*p+k*p*p]*Lambda[j+j*p+k*p*p]-
            Lambda[i+i*p+k*p*p]*Lambda[i+i*p+k*p*p];
      sl2 += (Lambda[j+j*p+k*p*p]-Lambda[i+i*p+k*p*p])*
             (Lambda[j+j*p+k*p*p]-Lambda[i+i*p+k*p*p]);
   
     }
   
     ASV += sl*sl*D_lm(F,p,q,i,j,0,0,B);

      for(k=0; k<K; k++){
        ASV += 2*sl*(Lambda[i+i*p+k*p*p]-Lambda[j+j*p+k*p*p])*
               D_lm(F,p,q,i,j,taus[k],0,B);
      } 

      ASV = ASV/(sl2*sl2);
  
      ASCOV -= ASV;
    
      for(k=0; k<K; k++){
        ASCOV += -((Lambda[i+i*p+k*p*p]-Lambda[j+j*p+k*p*p])*
                   D_lm(F,p,q,i,j,taus[k],0,B)-(Lambda[i+i*p+k*p*p]*
                  Lambda[i+i*p+k*p*p]-Lambda[j+j*p+k*p*p]*Lambda[i+i*p+k*p*p])*
                  D_lm(F,p,q,i,j,0,0,B))/sl2;
      }  

   
     result[0] = ASV;
     result[1] = ASCOV;
    

  }

  void ascov_deflji(double *F, double *Lambda, double *taus, int *nk, double *B, double *result)
  {
    int i=nk[0];
    int j=nk[1];
    int p=nk[2];
    int q=nk[3];
    int K=nk[4];
    int k;
    int l;
    int m;
 
    double ASV = 0;
    double ASCOV = 0;    
    double ujj = 0;
    double uij = 0;
    double sDk0 = 0;    

    for(k=0; k<K; k++){
      ASV += Lambda[j+j*p+k*p*p]*Lambda[j+j*p+k*p*p]*D_lm(F,p,q,i,j,taus[k],taus[k],B);
    }
 
    if(K>1){
    for(l=0; l<(K-1); l++){
      for(m=(l+1); m<K; m++){
          ASV += 2*Lambda[j+j*p+l*p*p]*Lambda[j+j*p+m*p*p]*D_lm(F,p,q,i,j,taus[l],taus[m],B);
        }
       }
      }      
     
     for(k=0; k<K; k++){
      ujj += Lambda[j+j*p+k*p*p]*Lambda[j+j*p+k*p*p];
      uij += Lambda[i+i*p+k*p*p]*Lambda[j+j*p+k*p*p];   
     }
   
     ASV += ujj*ujj*D_lm(F,p,q,i,j,0,0,B);

      for(k=0; k<K; k++){
        ASV -= 2*ujj*Lambda[j+j*p+k*p*p]*D_lm(F,p,q,i,j,taus[k],0,B);
      } 

      ASV = ASV/((ujj-uij)*(ujj-uij));
  
      ASCOV -= ASV;
    
      for(k=0; k<K; k++){
       sDk0 += Lambda[j+j*p+k*p*p]*D_lm(F,p,q,i,j,taus[k],0,B);
      
      }  

      ASCOV += (ujj*D_lm(F,p,q,i,j,0,0,B)-sDk0)/(ujj-uij);
   
     result[0] = ASV;
     result[1] = ASCOV;
  }

    void ascov_deflij(double *F, double *Lambda, double *taus, int *nk, double *B, double *result)
  {
    int i=nk[0];
    int j=nk[1];
    int p=nk[2];
    int q=nk[3];
    int K=nk[4];
    int k;
    int l;
    int m;
 
    double ASV = 0;
    double ASCOV = 0;    
    double uii = 0;
    double uij = 0;
    double sDk0 = 0;    

    for(k=0; k<K; k++){
      ASV += Lambda[i+i*p+k*p*p]*Lambda[i+i*p+k*p*p]*D_lm(F,p,q,i,j,taus[k],taus[k],B);
    }
 
    if(K>1){
    for(l=0; l<(K-1); l++){
      for(m=(l+1); m<K; m++){
          ASV += 2*Lambda[i+i*p+l*p*p]*Lambda[i+i*p+m*p*p]*D_lm(F,p,q,i,j,taus[l],taus[m],B);
        }
       }
      }      
     
     for(k=0; k<K; k++){
      uii += Lambda[i+i*p+k*p*p]*Lambda[i+i*p+k*p*p];
      uij += Lambda[i+i*p+k*p*p]*Lambda[j+j*p+k*p*p];   
     }
   
     ASV += uij*uij*D_lm(F,p,q,i,j,0,0,B);

      for(k=0; k<K; k++){
        ASV -= 2*uij*Lambda[i+i*p+k*p*p]*D_lm(F,p,q,i,j,taus[k],0,B);
      } 

      ASV = ASV/((uii-uij)*(uii-uij));
  
      ASCOV -= ASV;
    
      for(k=0; k<K; k++){
       sDk0 += Lambda[i+i*p+k*p*p]*D_lm(F,p,q,i,j,taus[k],0,B);
      }  

      ASCOV += (sDk0-uij*D_lm(F,p,q,i,j,0,0,B))/(uii-uij);
   
     result[0] = ASV;
     result[1] = ASCOV;
  }

  void D_lmij(double *F, int p, int q, int i, int j, int l, int m, double *B, double *result)
  { 
    result[0] = D_lm(F,p,q,i,j,l,m,B);
  }

 
}
