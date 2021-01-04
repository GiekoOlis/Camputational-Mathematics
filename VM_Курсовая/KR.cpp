#include <iostream>
#include <math.h>
using namespace std;


double a=3,b=4,c=2;
double EPS = 0.00000001;
double H=0.2,h=H/2;

double X_begin=0,X_end=1;
double Y_begin=1,YD_end=2.718281828459;
double YD_begin=-5;

double m1=YD_begin,m2=1.3,m3;

double *Y1,*Y2;
double shot[3];


double f1(double x, double y, double z) {
	return z;
}
double f2(double x, double y, double z) {
	return (a*exp(x)+b*z+c*y)/(a+b+c);
}
double RungeKuttFour(double h,int N,double *yy,double YD_begin){
	double X=X_begin, Y=Y_begin, Z=YD_begin;
	double K[4],L[4];
	for (int i=0;i<N;i++) {
		//cout<<"\t\t"<<X;
		//printf("\t\t%0.8f\t%0.8f\t\n",Y, Z);
		K[0]=h*f2(X, Y, Z);
		L[0]=h*f1(X, Y, Z);
		K[1]=h*f2(X+h/2, Y+L[0]/2, Z+K[0]/2);
		L[1]=h*f1(X+h/2, Y+L[0]/2, Z+K[0]/2);
		K[2]=h*f2(X+h/2, Y+L[1]/2, Z+K[1]/2);
		L[2]=h*f1(X+h/2, Y+L[1]/2, Z+K[1]/2);
		K[3]=h*f2(X+h, Y+L[2], Z+K[2]);
		L[3]=h*f1(X+h, Y+L[2], Z+K[2]);
		Z=Z+(K[0]+2*K[1]+2*K[2]+K[3])/6;		//y'
		yy[i]=Y+(L[0]+2*L[1]+2*L[2]+L[3])/6;	//y
		Y=yy[i];
		X+=h;
	}
	//cout<<"\t\t"<<X;
	//printf("\t\t%0.8f\t%0.8f\t\n",Y, Z);
	return Y;
}
double RK_DoubleRecount(double m) {
	int n1, n2;
	//cout << " \n\t\tX\t\tY\t\tY'" << endl;
	double result;
	do {
		n1=(X_end-X_begin)/H;
		n2=(X_end-X_begin)/h;
		Y1=new double[n1];
		Y2=new double[n2];
		//cout << "H = " << H << endl;
		RungeKuttFour(H, n1, Y1, m);
		//cout << "h = " << h << endl;
		result=RungeKuttFour(h, n2, Y2, m);
		//cout << endl;
		h*=0.5;
		H=2*h;
	} while(abs(Y2[n2-1]-Y1[n1-1])>EPS);
	h=0.1;
	H=2*h;
	return result;
}
void ShootingMethod() {
	cout<<"Первое приближение метода: "<<m1<<endl;
	shot[0]=RK_DoubleRecount(m1);
	cout<<"Второе приближение метода: "<<m2<<endl;
	shot[1]=RK_DoubleRecount(m2);
	if(abs(shot[0]-YD_end)<EPS){
		cout<<"Краевое условие для y' выполняется в точке x = "<<m1<<endl;
		printf("Решение краевой задачи: %0.8f\n\n", shot[0]);
		return;
	}
	else if(abs(shot[1]-YD_end)<EPS){
		cout<<"Краевое условие для y' выполняется в точке x = "<<m2<<endl;
		printf("Решение краевой задачи: %0.8f\n\n", shot[1]);
		return;
	}
	else {
		m3=m2+(((m2-m1)*(YD_end-shot[1]))/((shot[1]-shot[0])));
		cout<<"Новое приближение метода: "<<m3<<endl;
		shot[2]=RK_DoubleRecount(m3);
	}
	while (abs(shot[2]-YD_end)>=EPS) {
		m1=m2;
		m2=m3;
		shot[0]=shot[1];
		shot[1]=shot[2];
		m3=m2+(((m2-m1)*(YD_end-shot[1]))/((shot[1]-shot[0])));
		YD_begin=m3;
		cout<<"Новое приближение метода: "<<m3<<endl;
		shot[2]=RK_DoubleRecount(m3);
	}
	cout<<"Краевое условие для y' выполняется в точке x = "<<m3<<endl;
	printf("Решение краевой задачи: %0.8f\n\n", shot[2]);
	cout<<endl;
}
int main() {
	cout<<"\n\n\t\tРешение краевой задачи методом Рунге-Кутта IV порядка"<<endl;
	cout<<"Точность вычислений: "<<EPS<<endl;
	cout<<"Начальный шаг: "<<H<<endl;
	cout<<"Интервал: [ "<<X_begin<<" ; "<<X_end<<" ]"<<endl;
	cout<<"Выражение:"<<endl;
	cout<<"\t("<<a<<"exp(x)+"<<b<<"y'+"<<c<<"y)/"<<a+b+c<<endl;
	cout<<"Краевые условия:"<<endl;
	cout<<"\ty( "<<X_begin<<" ) = "<<Y_begin<<endl;
	printf("\ty'( %.0f ) = %.8f \n",X_end,YD_end);
	ShootingMethod();
	return 0;
}
