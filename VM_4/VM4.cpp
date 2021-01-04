#include <iostream>
#include <algorithm>
#include <math.h>
using namespace std;

#define n 2 //количество функций
#define m 2 //количество неизвестный
#define e 0.0001 //точность
#define h 1e-10

/*
#define prm1 0
#define prm2 -1
#define F1(x1,x2) sin(x1+x2)-1.6*x1
#define F2(x1,x2) pow(x1,2)+pow(x2,2)-1
*/

#define prm1 0.25 
#define prm2 0.75
#define F1(x1,x2) 0.1*pow(x1,2)+x1+0.2*pow(x2,2)-0.3
#define F2(x1,x2) 0.2*pow(x1,2)+x2-0.1*x1*x2-0.7


double func(int numb, double x1, double x2){
		switch(numb){
			case 1:
				return F1(x1,x2);
			case 2:
				return F2(x1,x2);
			default:
				return -1;
		}
}
double diff(int numb, int var, double x1, double x2){
	switch(var){
		case 1:
			return (func(numb,x1+h,x2)-func(numb,x1-h,x2))/(2.0*h);
		case 2:
			return (func(numb,x1,x2+h)-func(numb,x1,x2-h))/(2.0*h);
		default:
			return -1;
	}
		
}
double** YAkoby(double **J,double **X){
	// n - число строк и столбцов
    for(short int i=1;i<=n;i++){
		for(short int j=1;j<=n;j++){
			J[i-1][j-1]=diff(i,j,X[0][0],X[1][0]);
		}
	}
	
    return J;
}
double** InverseYAkoby(double **J){
    int leader_pos;
    double leader, temp;
	double **J_ = new double *[n];
	for (short int i=0;i<n;i++){
		J_[i] = new double[m];
	}     
	for(short int i=0;i<n;i++)
		for(short int j=0;j<n;j++)
			if (i==j) J_[i][j]=1;
    for(short int i=0;i<n;i++){
		leader=J[i][i];
		leader_pos=i;
		for(short int j=i;j<n;j++) {
			if ( fabs(J[j][i])>fabs(leader) ){
				leader=J[j][i];
				leader_pos=j;
			}    
		}
		if (fabs(leader)>h){
			for(short int j=0;j<n;j++){
				swap(J[i][j],J[leader_pos][j]);
				swap(J_[i][j],J_[leader_pos][j]);
			}
			for(short int j=0;j<n;j++){
				if (i!=j){
					temp=(J[j][i]/J[i][i]);
					for(short int k=0;k<n;k++){
						J[j][k]-=J[i][k]*temp;
						J_[j][k]-=J_[i][k]*temp;
					}
				}
			}
		}
	}  
	double t4;
	for(short int i=0;i<n;i++){
		t4=J[i][i];
		J[i][i]/=J[i][i]; 
		for(short int k=0;k<n;k++)
		J_[i][k]/=t4;
	}
return J_;
}
bool NormVector(double **X){
    double norm=0;
    double *Z = new double[n];
        for(short int j=0;j<n;j++){
            Z[j]=fabs(X[j][1]-X[j][0]);
            norm=max(norm,Z[j]);
        }
    if(norm>e) return true;
    return false;
}
int k=0;
double** Iteration(double **X){
	double *Z = new double[n];
	double sum=0;
	double **J = new double *[n];
		for (short int i=0;i<n;i++){
			J[i] = new double[m];
		}
	double F_x[n];
	for(short int j=0;j<n;j++) {
		X[j][0]=X[j][1];
		X[j][1]=0;
	}
	F_x[0]=func(1,X[0][0],X[1][0]);
	F_x[1]=func(2,X[0][0],X[1][0]);
	J=YAkoby(J,X);
	cout<<"i="<<k+1<<endl<<F_x[0]<<endl<<F_x[1]<<endl<<"Якоби: "<<endl;
	for(short int i=0;i<n;i++){
		for(short int j=0;j<n;j++){
			cout<<J[i][j]<<"\t";
		}
		cout<<endl;
	}
	J=InverseYAkoby(J);
	cout<<"Обратная: "<<endl;
	for(short int i=0;i<n;i++){
		for(short int j=0;j<n;j++){
			cout<<J[i][j]<<"\t";
		}
		cout<<endl;
	}
	for(short int i=0;i<n;i++){
		for(short int j=0;j<n;j++){
			sum+=J[i][j]*F_x[j];
		}
		Z[i]=sum; sum=0;
		X[i][1]=X[i][0]-Z[i];
	}
	cout<<"X="<<X[0][1]<<"\t"<<X[1][1]<<endl;
	cout<<"________"<<endl;
	return X;
}

double** MethodNewton(double **X){
	bool p; int max=25;
    X=Iteration(X);
    p=NormVector(X);
	k++;
    while(p==1 ){
		if(k>max) break;
        X=Iteration(X);
        p=NormVector(X);
		k++;
    }
	return X;
}

int main(){
	double **Result= new double *[n];
	for (short int i=0;i<n;i++) {
		Result[i] = new double[2];
	}  
	Result[0][1]=prm1; Result[1][1]=prm2;
	cout<<"Система уравнений:"<<endl<<"0.1*x1^2+x1+0.2*x2^2-0.3=0;"<<endl<<"0.2*x1^2+x2-0.1*x1*x2-0.7=0"<<endl;
	cout<<"Точность e=0.0001"<<endl;
	cout<<"Начальные приближения:\nX[1]="<<prm1<<endl<<"X[1]="<<prm2<<endl;
	Result=MethodNewton(Result);
	cout<<"Количество итераций: "<<k<<endl;
	cout<<"Ответ:"<<endl;
	for(short int j=0;j<m;j++){
		cout<<"X["<<j<<"]="<<Result[j][1]<<endl;
	}
}
