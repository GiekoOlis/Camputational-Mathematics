#include <iostream>
#include <algorithm>
#include <math.h>
#define ERROR -1
using namespace std;

double f(double x){
	return pow(x,0.333334); 
}

bool Check_Uniformity(double *X,double *F_x,int n){
	double h=X[1]-X[0];
	for(short int i=2;i<n;i++){
		if((X[i]-X[i-1])!=h) return false;
	}
	return true;
}
double Polinomial_Lagrange(double *X,double *F_x,double x,int n){
	double L_x,li_x=1;
	if(Check_Uniformity(X,F_x,n)){
		double h=X[1]-X[0];
		for (short int i=0;i<n;i++){
			for (short int j=0;j<n;j++){
				if (i!=j){
					li_x*=(x-X[0]-h*j)/h/(i-j);
				}
			}
			L_x+=F_x[i]*li_x;
			li_x=1;
		}
		return L_x;	
	}
	else{
		for(short int i=0;i<n;i++){
			for(short int j=0;j<n;j++){
				if(i!=j){
					li_x*=(x-X[j])/(X[i]-X[j]);
				}
			}
			L_x+=F_x[i]*li_x;
			li_x=1;
		}
		return L_x;
	}
}
double Aitkens_Scheme(double *X, double *F_x, double x, int n){
	double Z[n];
	for(short int i=0;i<n;i++){
		Z[i]=F_x[i];
	}
	for(short int j=0;j<n;j++){
		for(short int i=j+1;i<n;i++){
			Z[i]=((x-X[i])*Z[j]-(x-X[j])*Z[i])/(X[i]-X[j]);
		}
	}
	return Z[n-1];
}

double Polinomial_Newton(double *X, double *F_x, double x, int n){
	double *k=new double[n];
	k[0]=F_x[0];
	for (int j=1; j<n; j++)
		for (int i=0; i<n-j; i++){
			F_x[i]=(F_x[i+1]-F_x[i])/(X[i+j]-X[i]);
			k[j]=F_x[0];
		}
	double S=k[0], p=1;
	for (int i=1; i<n; i++){
		p*=(x-X[i-1]);
		S+=k[i]*p;
	}
	return S;
}


int main(){

	int n=5; double x=1.5; double Res;
	double *F_x= new double [n];
	double X[n]={1,1.3,1.6,1.8,2.1};
	
	
	cout<<"F(x)=pow(x,0.333334)"<<endl;
	cout<<"x\tF(x)"<<endl;
	for(short int i=0;i<n;i++){
		F_x[i]=f(X[i]);
		cout<<X[i]<<"\t"<<F_x[i]<<endl;
	}
	cout<<"X="<<x<<endl<<"F(X)=F("<<x<<")="<<f(x)<<endl<<endl;
	Res=Polinomial_Lagrange(X,F_x,x,n);
	cout<<"Polinomial Lagrange"<<endl;
	cout<<"Result: "<<Res<<endl<<endl;
	Res=Aitkens_Scheme(X,F_x,x,n);
	cout<<"Aitkens_Scheme"<<endl;
	cout<<"Result: "<<Res<<endl<<endl;
	Res=Polinomial_Newton(X,F_x,x,n);
	cout<<"Polinomial Newton"<<endl;
	cout<<"Result: "<<Res<<endl;
	return 0;
}



