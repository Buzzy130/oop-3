#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "function.h"
#include "empirical_main.h"
#include "mixture_main.h"
#include "huber_main.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <string>
using namespace std;

void next()
{
	char option;

	cout << "Распределение: " << endl;
	cout << "1. Основное распределение" << endl;
	cout << "2. Смесь распределений" << endl;
	cout << "3. Эмпирическое распределение" << endl;
	cout << "4. Выход" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		Huber_D();
		break;

	case '2':
		mixture_D();
		break;

	case '3':
		empirical_D();
		break;

	case '4':
		cout << "Завершение работы" << endl << endl;
		exit(0);

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}


int main(int argc, char** argv)
{
	setlocale(LC_ALL, "ru");
	int n = 1000;
	//int result = Catch::Session().run(argc, argv);
	//return result;
	next();
}

