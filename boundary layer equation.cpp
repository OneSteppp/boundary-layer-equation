#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;
int main()
{
	double Max_Eta = 8;//Eta范围
	int N = 1000;//节点个数
	double h = Max_Eta / N;//步长
	double *f1, *f2, *f3, *f4, *f5, *f6, *Eta;
	double K1, K2, K3, K4, L1, L2, L3, L4, M1, M2, M3, M4, N1, N2, N3, N4, P1, P2, P3, P4, Q1, Q2, Q3, Q4;//龙格库塔参数
	double S = 1.1, S_new = 1, S_old=0, m = 0;//m=-0.05,0,0.2,1
	int i;

	f1 = new double[N + 1];
	f2 = new double[N + 1];
	f3 = new double[N + 1];
	f4 = new double[N + 1];
	f5 = new double[N + 1];
	f6 = new double[N + 1];
	Eta = new double[N + 1];//分配空间

	f1[0] = 0, f2[0] = 0, f3[0] = S, f4[0] = 0, f5[0] = 0, f6[0] = 1;//定义初值

	for (i = 0; i <= N; i++)
	{
		Eta[i] = i*h;
	}//横坐标

	do
	{
		f3[0] = S;
		for (i = 0; i<N; i++)
		{
			K1 = f2[i];
			L1 = f3[i];
			M1 = -0.5*(m + 1)*f1[i] * f3[i] - m*(1 - f2[i] * f2[i]);
			N1 = f5[i];
			P1 = f6[i];
			Q1 = -0.5*(m + 1)*(f3[i] * f4[i] + f1[i] * f6[i]) + 2 * m*f2[i] * f5[i];
			//
			//
			K2 = f2[i] + 0.5*h*L1;
			L2 = f3[i] + 0.5*h*M1;
			M2 = -0.5*(m + 1)*(f1[i] + 0.5*h*K1)*(f3[i] + 0.5*h*M1) - m*(1 - (f2[i] + 0.5*h*L1)*(f2[i] + 0.5*h*L1));
			N2 = f5[i] + 0.5*h*P1;
			P2 = f6[i] + 0.5*h*Q1;
			Q2 = -0.5*(m + 1)*((f3[i] + 0.5*h*M1)*(f4[i] + 0.5*h*N1) + (f1[i] + 0.5*h*K1)*(f6[i] + 0.5*h*Q1)) + 2 * m*(f2[i] + 0.5*h*L1)*(f5[i] + 0.5*h*P1);
			//
			//
			K3 = f2[i] + 0.5*h*L2;
			L3 = f3[i] + 0.5*h*M2;
			M3 = -0.5*(m + 1)*(f1[i] + 0.5*h*K2)*(f3[i] + 0.5*h*M2) - m*(1 - (f2[i] + 0.5*h*L2)*(f2[i] + 0.5*h*L2));
			N3 = f5[i] + 0.5*h*P2;
			P3 = f6[i] + 0.5*h*Q2;
			Q3 = -0.5*(m + 1)*((f3[i] + 0.5*h*M2)*(f4[i] + 0.5*h*N2) + (f1[i] + 0.5*h*K2)*(f6[i] + 0.5*h*Q2)) + 2 * m*(f2[i] + 0.5*h*L2)*(f5[i] + 0.5*h*P2);
			//
			//
			K4 = f2[i] + h*L3;
			L4 = f3[i] + h*M3;
			M4 = -0.5*(m + 1)*(f1[i] + h*K3)*(f3[i] + h*M3) - m*(1 - (f2[i] + h*L3)*(f2[i] + h*L3));
			N4 = f5[i] + h*P3;
			P4 = f6[i] + h*Q3;
			Q4 = -0.5*(m + 1)*((f3[i] + h*M3 )*(f4[i] + h*N3) + (f1[i] + h*K3)*(f6[i] + h*Q3)) + 2 * m*(f2[i] + h*L3)*(f5[i] + h*P3);
			//
			f1[i + 1] = f1[i] + h*(K1 + 2*K2 + 2*K3 + K4) / 6.0;
			f2[i + 1] = f2[i] + h*(L1 + 2*L2 + 2*L3 + L4) / 6.0;
			f3[i + 1] = f3[i] + h*(M1 + 2*M2 + 2*M3 + M4) / 6.0;
			f4[i + 1] = f4[i] + h*(N1 + 2*N2 + 2*N3 + N4) / 6.0;
			f5[i + 1] = f5[i] + h*(P1 + 2*P2 + 2*P3 + P4) / 6.0;
			f6[i + 1] = f6[i] + h*(Q1 + 2*Q2 + 2*Q3 + Q4) / 6.0;
			//龙格-库塔法
		}
		/*
		     f1'=f2, f2'=f3, f3'=-0.5*(m+1)*f1*f3-m*(1-f2*f2)
		*/
		/*
		     f4=df1/dS, f5=df2/dS, f6=df3/dS
			 f4'=f5, f5'=f6, f6'=-0.5*(m+1)(f1*f6+f3*f4)+2*m*f2*f5
	    */
		S_old = S;
		S_new = S + (1 - f2[N]) / f5[N];
		S = S + (1 - f2[N]) / f5[N];
		/*
		f5=df2/dS; (f2[N](应有值)-f2[N])/(S_应有值-S）=f5[N];
		S_应有值=S+(1-f2[N])/f5[N]
		*/
		cout << S << endl;
	} while (fabs(S_new - S_old) > 0.001);
	//
	//输出
	ofstream ofile("result_1.txt");
	for (int i = 0; i<= N; i++)
	{
		ofile << Eta[i] << '\t' << f2[i] << endl;;
	}
	ofile.close();
	delete[] f1;
	f1 = NULL;
	delete[] f2;
	f2 = NULL;
	delete[] f3;
	f3 = NULL;
	delete[] f4;
	f4 = NULL;
	delete[] f5;
	f5 = NULL;
	delete[] f6;
	f6 = NULL;
	delete[] Eta;
	Eta = NULL;
	return 0;
}
