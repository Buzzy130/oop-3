#include "huber_main.h"

void end_huber(HuberD* HD)
{
	ofstream x_selection;
	ofstream y_selection;
	ofstream params;
	vector<pair<double, double>> table;

	int n;
	double x = 0;

	cout << " Вывести параметры распределения на экран" << endl;
	cout << "v = " << HD->get_v() << endl;
	cout << "K = " << HD->get_k() << endl;
	cout << "Scale = " << HD->get_scale() << endl;
	cout << "Shift = " << HD->get_shift() << endl;
	cout << "M[X] = " << HD->Mksi_huber() << endl;
	cout << "D[X] = " << HD->Dksi_huber() << endl;
	cout << "Asymmetry[X] = " << HD->asymmetry_huber() << endl;
	cout << "Kurtosis[X] = " << HD->kurtosis_huber() << endl;
	cout << "P = " << HD->P() << endl;
	cout << "Density(0) = " << HD->Huber(0) << endl << endl << endl;

	cout << "Вычислить значение плотности распределения в произвольной точке" << endl;
	cout << "Введите x: ";
	cin >> x;
	cout << endl << "Density(" << x << ") = " << HD->Huber(x) << endl << endl;

	x_selection.open("selection1.txt");
	y_selection.open("selection2.txt");

	cout << "Получить выборку для анализа" << endl;
	cout << "Введите размерность выборки n: ";
	cin >> n;
	table = HD->generate_pair(n);
	for (const pair<double, double>& pr : table)
	{
		x_selection << pr.first << endl;
		y_selection << pr.second << endl;
	}

	x_selection.close();
	x_selection.close();
	cout << " Вывести параметры распределения в файл" << endl;
	cout << endl << "Значения выборки записаны в файлы selection1.txt и selection2.txt" << endl << endl;

		HD->save_file(params);
		cout << endl << endl;	
}

void Huber_D()
{
	char option;
	double v, scale, shift;
	HuberD* HD;
	ifstream file;

	cout << "Выберите опцию: " << endl;
	cout << "1. Стандартное распределение Хьюбера" << endl;
	cout << "2. Ввести произвольные параметры с клавиатуры" << endl;
	cout << "3. Считать параметры из файла" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "Введите параметр v: ";
		cin >> v;
		cout << endl << endl;
		try
		{
			HD = new HuberD(v);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}

		end_huber(HD);
		break;

	case '2':
		cout << "Введите v scale shift: ";
		cin >> v >> scale >> shift;
		cout << endl << endl;
		try
		{
			HD = new HuberD(v, scale, shift);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}

		end_huber(HD);
		break;

	case '3':
		try
		{
			HD = new HuberD(file);
		}
		catch (exception e)
		{
			cerr << endl << e.what() << endl << endl;
			break;
		}

		end_huber(HD);
		break;

	case '4':
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}