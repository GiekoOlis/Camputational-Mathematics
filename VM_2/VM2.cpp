#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
using namespace std;
//Функция для проверки сходимости матрицы по трем основным условиям
double ConvergenceMatrix(double **C,int n){
    double Norm1,Norm2,Infinity;
    double sum=0;
    for(short int j=0; j<n; j++){
        for (short int k=0; k<n; k++) sum+=abs(C[j][k]);
        Infinity = max(Infinity,sum);
        sum=0;
    }
    if(Infinity<1) return Infinity;
    else{
        for(short int j=0; j<n; j++){
            for (short int k=0; k<n; k++) sum+=abs(C[k][j]);
            Norm1 = max(Norm1,sum);
            sum=0;
            }
        if(Norm1<1) return Norm1;
        else{
            for(short int j=0; j<n; j++){
                for (short int k=0; k<n; k++) sum+=pow(C[k][j],2);
            }
            Norm2=sqrt(sum);
            if(Norm2<1) return Norm2;
            else return 0;
        }
    }
}
//Здесь проверяется, достигли ли мы заявленной точности
bool NormVector(double **X,int n,double e){
    double norm=0;
    double *Z = new double[n];
        for(short int j=0;j<n;j++){
            Z[j]=fabs(X[j][1]-X[j][0]);
            norm=max(norm,Z[j]);
        }
    if(norm>=e) return true;
    return false;
}
//Одна итерация метода итераций
//Включающая в себя сдвиг текущего и прошлого ответа,
//Поиск произведения A*X(i-1) для итерационной формулы
//Получение новых текущих результатов Х(i)
double** Iteration(double **C, double *B, double **X,int n){
     double *Z = new double[n];
     double sum=0;
     for(short int j=0;j<n;j++) {
            X[j][0]=X[j][1];
            X[j][1]=0;
         }
     for(short int i=0;i<n;i++){
         for(short int j=0;j<n;j++){
                sum+=C[i][j]*X[j][0];
            }
            Z[i]=sum; sum=0;
             X[i][1]=B[i]+Z[i];
        }
        return X;
}
//Поиск необходимого количества итераций 
//Априорная оценка 
int NumberIteration(double *B,double e,double q,int z){
    int n;double l;
    double norm=0;
    for (short int j=0; j<z; j++){
                norm = max(norm,abs(B[j]));
    }
    l=(log((e*(1-q))/(norm)))/(log(q));
    n=round(l);
    return n;
}
/*SimpleCalculate и HardCalculate - метод простых итераций*/
//Метод решения, основанный на заранее подсчитанном количестве необходимых итераций
double** SimpleCalculate(double **C,double *B,double **X,int n,double norm,double e){
    int k;
    k=NumberIteration(B,e,norm,n);
    cout<<"norm: "<<norm<<"\nnumber iteration: "<<k<<endl;
    for(short int i=0;i<k;i++){
    X=Iteration(C,B,X,n);
    }
    return X;
}
//Метод решения, основанный на проверке соответствия заявленной точности после каждой итерации
double** HardCalculate(double **C,double *B,double **X,int n,double e){
    bool p;
    X=Iteration(C,B,X,n);
    p=NormVector(X,n,e);
    while(p==1){
        X=Iteration(C,B,X,n);
        p=NormVector(X,n,e);
    }
    return X;
}
//Получение новых текущих результатов Х(i) для Зейделя
double** IterationSeidel(double **C1,double **C2, double *B, double **X,int n){
	double *Z1 = new double[n];
	double *Z2 = new double[n];
	double sum1=0,sum2=0;
	for(short int j=0;j<n;j++) {
            X[j][0]=X[j][1];
            X[j][1]=0;
         }
     for(short int i=0;i<n;i++){
         for(short int j=0;j<n;j++){
                sum1+=C1[i][j]*X[j][1];
				sum2+=C2[i][j]*X[j][0];
            }
            Z1[i]=sum1; sum1=0;
			Z2[i]=sum2; sum2=0;
			X[i][1]=B[i]+Z1[i]+Z2[i];
        }
	return X;
}
int k=0;
//Метод Зейделя 
double** Seidel(double **C,double *B,double **X, int n,double e){
	double **C1 = new double *[n];
	double **C2 = new double *[n];
	for (int i = 0; i < n; i++) {
		C1[i] = new double[n];
		C2[i] = new double[n];
	} 
	for(short int i=0;i<n;i++){
		for(short int j=0;j<n;j++){
			if(j>i){
				C2[i][j]=C[i][j];
				C1[i][j]=0;
			}
			else{
				C1[i][j]=C[i][j];
				C2[i][j]=0;
			}
		}
	}
	bool p;
	X=IterationSeidel(C1,C2,B,X,n);
    p=NormVector(X,n,e);k=1;
    while(p==1){
    X=IterationSeidel(C1,C2,B,X,n);    
	p=NormVector(X,n,e);
    k++;
    }
return X;
}




int main(){
int coef,numb;
double e;
ifstream in("testsVM2.txt");
if (in.is_open()) {
    in>>coef; 
    for (short int i=0; i<coef; i++){
        cout<<"---------------------------------"<<endl;
        cout<<"MATRIX "<<i+1<<endl;
        cout<<"---------------------------------"<<endl;
        in>>numb; in>>e;
        cout<<"e = "<<e<<endl;
        double *B = new double [numb];
        double **X = new double *[numb];
        double **C = new double *[numb];
        for (int i = 0; i < numb; i++) {
            C[i] = new double[numb];
            X[i] = new double[2];
        }            
        for (short int j=0; j<numb; j++){
            for (short int k=0; k<numb; k++){
                in>>C[j][k];
                cout<<C[j][k]<<"\t";
            }
            in>>B[j];
            cout<<"| "<<B[j]<<endl;
        }          
        for (short int j=0; j<numb; j++){
            double p=C[j][j];
            for (short int k=0; k<numb; k++){
                C[j][k]=-C[j][k]/p;
            }
            C[j][j]=0;
            B[j]=B[j]/p;
        }
    double norm=ConvergenceMatrix(C,numb);
    if(norm!=0)
        
    X=Seidel(C,B,X,numb,e);
    //X=HardCalculate(C,B,X,numb,e);
    //X=SimpleCalculate(C,B,X,numb,norm,e);
    else cout<<"Matrix not convergence"<<endl;
    cout<<"result"<<endl;
        for(short int j=0; j<numb; j++){
            cout<<"X["<<j+1<<"] = "<<X[j][1]<<endl; cout<<k<<endl;
        }
    }
    in.close();
}
else {
cout << "Not true way to file.."<<endl;
}
return 0;
}
