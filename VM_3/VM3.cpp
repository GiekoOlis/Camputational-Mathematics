#include <iostream>
#include <algorithm>
#include <math.h>
using namespace std;


double f(double x){
	return pow(x,3)-3; //x=1,44225
	//return pow(x,2)+4*x; //x=-4;x=0
	//return pow(x,3)+4*x+5; //x=-1
}
double diff(double x){
		const double h=1e-10;
	return (f(x+h)-f(x-h))/(2.0*h);
}
double Bisection(double a, double b,double e,int *k){
	double f_x,x,tmp;
	(*k)=1;
	x=(b+a)/2;
	tmp=f(a);  
	f_x=f(x);
	if(tmp*f_x<0) b=x;
	else a=x;
	while(fabs(f_x)>=e){
		(*k)++;
		x=a+(b-a)/2;
		tmp=f(a);  
		f_x=f(x);
		if(tmp*f_x<0) b=x;
		else a=x;
	}
	return x;
}
double Chord(double a, double b,double e,int *k){
	double f_x,x,tmp;
	x=(a*f(b)-b*f(a))/(f(b)-f(a));
	tmp=f(a);  
	f_x=f(x);
	(*k)=1;
	if(tmp*f_x<0) b=x;
	else a=x;
	while(fabs(f(x))>=e){
		(*k)++;
		x=(a*f(b)-b*f(a))/(f(b)-f(a));
		if((f(a)*f(x))<0) b=x;
		else a=x;
	}
	return x; 
}
double Newton(double a, double e, int *k){
	double x[2];
	(*k)=1;
	x[0]=a;
	x[1]=x[0]-f(x[0])/diff(x[0]);
	if ((fabs(x[1]-x[0])<=e)) return x[1];
	while(fabs(x[1]-x[0])>e){
		x[0]=x[1];
		x[1]=x[0]-f(x[0])/diff(x[0]);
		(*k)++;
	}
	return x[1]; 
}

int main(){
	double a=1.4,b=1.5,X;
	double e=0.00001;int k=0;
	cout<<"f(x)=x^3-3=0\nx=3^(1/3)\nx=1.44225"<<endl;
	cout<<"Интервал поиска: [1.4;1.5]"<<endl;
	cout<<"E = "<<e<<endl;
	cout<<endl;
	X=Bisection(a,b,e,&k);
	cout<<"Bisection: "<<endl;
	cout<<"X = "<<X<<endl;
	cout<<"Количество итераций: "<<k<<endl;
	cout<<endl;
	cout<<"Chord: "<<endl;
	X=Chord(a,b,e,&k);
	cout<<"X = "<<X<<endl;
	cout<<"Количество итераций: "<<k<<endl;
	cout<<endl;
	cout<<"Newton: "<<endl;
	X=Newton(a,e,&k);
	cout<<"X = "<<X<<endl;
	cout<<"Количество итераций: "<<k<<endl;
	return 0;   
}
