#include <iostream>
#include <algorithm>
#include <math.h>
#define ERROR -1
#define EPS 0.0000001
using namespace std;
double h=1e-10;
int pos;
double f(double x){
	return pow(x,0.333334); 
}
double** Transform(double **C,int n,int k){
	double max; 
	short int in=k;
	max=abs(C[k][k]);
	for(short int i=k+1;i<n;i++){
		if(abs(C[i][k])>max){
		max=abs(C[i][k]);
		in=i;
		}
	}
    if(max<EPS){
		cout << "null column"<<endl;
		return NULL;
	}
	for(short int j=0;j<=n;j++) {
		double temp=C[k][j];
		C[k][j]=C[in][j];
		C[in][j]=temp;
	}
	return C;
}
double* Gauss(double **C, double *X, int n){
	int k=0;
	while(k<n){
		C=Transform(C,n,k);
		if(C==NULL) return NULL;
		else{
			for(int i=k;i<n;i++){
				double temp=C[i][k];
				if (abs(temp) < EPS) continue;
				for(short int j=0;j<=n;j++) 
					C[i][j]/=temp;
				if(i==k)  continue;
				for(short int j=0;j<=n;j++)
					C[i][j]-=C[k][j];
			}
		}
		k++;	
	}
	for(k=n - 1;k>=0;k--){
		X[k]=C[k][n];
		for(int i=0;i<k;i++)
			C[i][n]-=C[i][k]*X[k];
	}
	return X;
}
int FindX(double *X,double x,int n){
	for(short int i=0;i<n;i++)
		if((x>=X[i])&&(x<=X[i+1])) return i+1;
	return ERROR;
}
double Cubic_Spline_Interpolation(double *X, double *A, double x,int n){
	double S;
    double *H= new double [n];
	for(short int i=1;i<n;i++) H[i]=X[i]-X[i-1];
	double *M = new double [n];
	double **C = new double *[n-2];
	for (short int i=0;i<n-2;i++) {
		C[i] = new double[n-1];
	}
	for (short int j=1; j<n-1; j++){
		for (short int k=1; k<n-1; k++){
			if(j==k)			C[j-1][k-1]=(H[j]+H[j+1])/3;
			else if(k==j+1)		C[j-1][k-1]=(H[j+1])/6;
			else if(k==j-1)		C[j-1][k-1]=H[j]/6;
			else 				C[j-1][k-1]=0;
		}
	}
	
	for (short int i=0; i<n-2; i++){
		C[i][n-2]=(A[i+2]-A[i+1])/H[i+2]-(A[i+1]-A[i])/H[i+1];
	}
	M=Gauss(C,M,n-2);
	for (short int k=n-2; k>0; k--)		M[k]=M[k-1];
	M[0]=0;M[n-1]=0;
	pos=FindX(X,x,n);
	S=0;
	S+=M[pos-1]*pow(X[pos]-x,3)/(6*H[pos]);
	S+=M[pos]*pow(x-X[pos-1],3)/(6*H[pos]);
	S+=(A[pos-1]-M[pos-1]*pow(H[pos],2)/6)*(X[pos]-x)/H[pos];
	S+=(A[pos]-M[pos]*pow(H[pos],2)/6)*(x-X[pos-1])/H[pos];
    return S;
}

int main(){
    int n=5; double x=2;
	double *F_x= new double [n];
	double X[n]={1,1.3,1.6,1.8,2.1};
	cout<<"F(x)=pow(x,0.333334)"<<endl;
	cout<<"x\tF(x)"<<endl;
	for(short int i=0;i<n;i++){
		F_x[i]=f(X[i]);
		cout<<X[i]<<"\t"<<F_x[i]<<endl;
	}
	cout<<"X="<<x<<endl<<"F(X)=F("<<x<<")="<<f(x)<<endl;
	double S=Cubic_Spline_Interpolation(X,F_x,x,n);
	cout<<"X="<<x<<" находится внутри отрезка ["<<X[pos-1]<<";"<<X[pos]<<"]"<<endl<<endl;
	cout<<"Cubic Spline Interpolation:"<<endl;
	cout<<"S["<<pos<<"][X]=S["<<pos<<"]["<<x<<"]="<<S<<endl;
	return 0;
}



