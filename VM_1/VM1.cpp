#include <iostream>
#include <fstream>
#include <math.h>
#define EPS 0.00000001
using namespace std;


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
int main(){
int coef,numb;//coef - number matrix
ifstream in("testsVM1.txt"); //classwork
if (in.is_open()) {
	in>>coef;
	for (short int i=0; i<coef; i++){
		cout<<"---------------------------------"<<endl;
		cout<<"MATRIX "<<i+1<<endl;
		cout<<"---------------------------------"<<endl;
		in>>numb;
		double *X = new double [numb];
		double **C = new double *[numb];
		for (int i = 0; i < numb; i++) {
			C[i] = new double[numb+1];
		}       
		for (short int j=0; j<numb; j++){
			for (short int k=0; k<numb; k++){
				in>>C[j][k];
				cout<<C[j][k]<<"\t";
			}
			in>>C[j][numb];
			cout<<"| "<<C[j][numb]<<endl;
		}
		X=Gauss(C,X,numb);
		if(X==NULL) cout<<"matrix can not be resolved"<<endl;
		else{
			for(short int j=0; j<numb; j++){
				cout<<"X["<<j+1<<"] = "<<X[j]<<endl;
			}
		}
	}
	in.close();
}
else {
cout << "Not true way to file.."<<endl;
}
return 0;
}
