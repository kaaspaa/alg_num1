#define _CRT_SECURE_NO_DEPRECATE
#define _USE_MATH_DEFINES
#include<stdio.h>
#include<iostream>
#include<cmath>
#include<fstream>
#define MAX 150
//funkcja : sin(x) * ln(1+x);
using namespace std;

double pwr(double x, long long p) {
	// x - zmienna, p- potega
	double answ = 1.0;
	if (p<0) {
		x = 1.0 / x;
		p = p * -1;
	}
	if (p == 0)
		answ = 1.0;
	if (p>0) {
		long long i;
		for (i = 1;i <= p;i++) {
			answ = answ * x;
		}
	}
	return answ;
}

double factorial(double zmienna) {
	// silnia
	double answ = 1;
	if (zmienna <0)
		answ = 0;
	else if (zmienna == 0)
		answ = 1;
	else {
		for (double i = 2.0;i <= zmienna;i++) {
			answ = answ * i;
		}
	}
	return answ;
}
//2(n-1) = 2n -1
double sin_front(double x, int n){
		double result = 0.0;
		for (int i = 0; i < n; i++)
		{
			result += (pwr(-1, i) * pwr(x, 2 * i + 1)) / factorial(2 * i + 1);
		}
		return result;
	}

double sin_back(double x, int n) {
	double result = 0.0;
	double tab[n];
	for (int i = 0; i < n; i++)
	{
		tab[i] = (pwr(-1, i) * pwr(x, 2 * i + 1)) / factorial(2 * i + 1);
	}
	for (int i = n - 1; i >= 0; i--)
	{
		result += tab[i];
	}
	return result;
}


static double sin_front_on_prev(double x, int n)
{
	double result = x;
	double tab[n];
	tab[0] = x;
	for (int i = 1; i < n; i++)
	{
		tab[i] = (tab[i - 1] * (-1)*(pwr(x, 2))) / (2 * i*(2 * i + 1));
		result += tab[i];
	}

	return result;
}
static double sin_back_on_prev(double x, int n)
{
	double tab[n];
	tab[0] = x;
	for (int i = 1; i < n; i++)
	{
		tab[i] = (tab[i - 1] * (-1)*(pwr(x, 2))) / (2 * i*(2 * i + 1));
	}
	double result = 0.0;
	for (int i = n - 1; i >= 0; i--)
	{
		result += tab[i];
	}
	return result;
}

double cos_front(double x, int n) {
	double result=0.0;
	for (int i = 0;i < n;i++) {
		result = result + (pwr(-1, i) * pwr(x, 2 * i)) / factorial(2 * i);
	}
	return result;
}

double cos_back(double x, int n) {
	double result = 0.0;
	double tab[n];
	for (int l = 0;l < n;l++) {
		tab[l] = (pwr(-1, l) * pwr(x, 2 * l)) / factorial(2 * l);
	}
	for (int i = n-1;i >= 0;i--) {
		result = result + tab[i];
	}
	return result;
}

static double cos_front_on_prev(double x, int n) {
	double result = 0.0;
	double tab[n];
	tab[0] = 1;
	for (int i = 1; i < n; i++)
		tab[i] = (tab[i - 1] * (-1)*(pwr(x, 2))) / ((2 * i)*(2 * i - 1));
	for (int l = 0;l < n;l++)
		result = result + tab[l];
	return result;
}

static double cos_back_on_prev(double x, int n) {
	double result = 0.0;
	double tab[n];
	tab[0] = 1;
	for (int i = 1; i < n; i++)
		tab[i] = (tab[i - 1] * (-1)*(pwr(x, 2))) / ((2 * i)*(2 * i - 1));
	for (int l = n-1;l >= 0;l--)
		result = result + tab[l];
	return result;
}

int main() {
	//x - argumenty, n - zasieg
	double x;
	int n = 5;
	//
	double start_val = M_PI * (-1);
	double end_val = M_PI;
	double co_ile = (end_val - start_val) / 1000000;
	//lim - 1000 zeby bylo wiadomo kiedy przerwac liczyc i dac do wykresu
	int lim = 100;
	//pomocnicze
	int i,l,o;
	//plik/pomocnicze/odpowiedz_z_cmath
	FILE * f;
	f = fopen("Wyniki.csv", "w+");
	//tablica bledow i jej zerowanie
	long double blad[4][lim];
	double przedzial = 1000000/lim;
	for(i=0;i<lim;i++){
		for(l=0;l<4;l++)
			blad[l][i] = 0.0;
	}
	i=0;
	l=0;
	o=0;//pomocnicza do przedzialu
	long double lib_answ;//odpowiedz z biblioteki
	for(x=start_val;x<end_val;x=x+co_ile){
		lib_answ = sin(cos(x));
		blad[l][i] += fabs(lib_answ - sin_front(cos_front(x,n),n));
		l++;
		blad[l][i] += fabs(lib_answ - sin_front_on_prev(cos_front_on_prev(x,n),n));
		l++;
		blad[l][i] += fabs(lib_answ - sin_back(cos_back(x,n),n));
		l++;
		blad[l][i] += fabs(lib_answ - sin_back_on_prev(cos_back_on_prev(x,n),n));
		l=0;
		o++;
		if(o == przedzial){
			i++;
			o=0;
		}
	}
	double wyn[4];
	for(l=0;l<4;l++)
		wyn[l] = 0.0;//zerowanie tablicy wynikow	

	for(l=0;l<4;l++){
		for(i=0;i<lim;i++){
			blad[l][i] = blad[l][i] / przedzial;//usrednianie
			wyn[l] = wyn[l] + blad[l][i];//dodawanie usrednionych wynikow do tablicy
			fprintf(f,"%.20Lf\n", blad[l][i]);//wypisywanie
		}
	fprintf(f,"\n");
	}
	//zamkniecie pliku
	fclose(f);
	return 0;
}
