#include "empirical_main.h"


void end_empirical(Empirical* ED)
{
	ofstream file;
	double x = 0;

	cout << "Вывести параметры эмпирического распределения на экран" << endl;
	cout << "n = " << ED->get_n() << endl;
	cout << "K = " << ED->get_k() << endl;
	cout << "M[X] = " << ED->Mn() << endl;
	cout << "D[X] = " << ED->Dn() << endl;
	cout << "As[X] = " << ED->asymmetry_empirical() << endl;
	cout << "Kurt[X] = " << ED->kurtosis_empirical() << endl;
	cout << "f(0) = " << ED->H_Empirical(x) << endl << endl << endl;

	cout << "Вычислить значение плотности смеси распределений в произвольной точке" << endl;
	cout << "Введите x: ";
	cin >> x;
	cout << endl << "f(" << x << ") = " << ED->H_Empirical(x) << endl << endl;

	cout << "3. Вывести выборку элементов и параметры эмпирического распределения в файл" << endl;
	ED->save_to_file(file);
	cout << endl << endl;

}

void empirical_huber()
{
	Empirical* ED;
	char option;
	double v, scale, shift;
	int n, k;
	HuberD* HD;
	ifstream file;

	cout << "Выберите опцию: " << endl;
	cout << "1. Стандартное распределение Хьюбера" << endl;
	cout << "2. Ввести произвольные параметры с клавиатуры" << endl;
	cout << "3. Считать параметры из файла" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl << endl;

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
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		ED = new Empirical(HD, n, k);
		end_empirical(ED);

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
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;
		ED = new Empirical(HD, n, k);
		end_empirical(ED);

		break;

	case '3':
		try
		{
			HD = new HuberD(file);
		}
		catch (exception e)
		{
			cerr << endl << e.what() << endl;
			break;
		}

		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;

		ED = new Empirical(HD, n, k);
		end_empirical(ED);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}

void empirical_mixture()
{
	Empirical* EM;
	char option;
	double p, v1, scale1, shift1, v2, scale2, shift2;
	int n, k;
	Mixture* MD;
	HuberD* HD1;
	HuberD* HD2;
	ifstream file;

	cout << "Выберите опцию: " << endl;
	cout << "1. Стандартные параметры смеси распределений" << endl;
	cout << "2. Ввод с клавиатуры" << endl;
	cout << "3. Ввод из файла" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		MD = new Mixture();

		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MD, n, k);

		end_empirical(EM);
		break;

	case '2':
		cout << "Введите p (v1 scale1 shift1) (v2 scale2 shift2) (через пробел или через Enter): ";
		cin >> p >> v1 >> scale1 >> shift1 >> v2 >> scale2 >> shift2;
		try
		{
			HD1 = new HuberD(v1, scale1, shift1);
			HD2 = new HuberD(v2, scale2, shift2);
			MD = new Mixture(HD1, HD2, p);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MD, n, k);
		end_empirical(EM);
		break;

	case '3':
		try
		{
			MD = new Mixture(file);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}

		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MD, n, k);
		end_empirical(EM);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}

void empirical_empirical()
{
	Empirical* EM;
	char option;
	ifstream file;
	string filename;
	int n, k;

	cout << "1. Ввести параметры распределения из файла" << endl;
	cout << "2. Назад" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		cout << "Введите имя файла с параметрами распределения: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("Ошибка: не удалось открыть файл");

		while (!file.eof())
			file >> n >> k;

		file.close();

		try
		{
			EM = new Empirical(n, k);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		cout << endl << endl;

		end_empirical(EM);
		break;

	case '2':
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;

	}
}

void empirical_selection()
{
	Empirical* EM;
	char option;
	ifstream file;
	string filename;
	double x;
	vector<double> x_s;

	cout << "1. Ввести выборку элементов из файла" << endl;
	cout << "2. Назад" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "Введите имя файла с параметрами распределения: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("Ошибка: не удалось открыть файл");

		while (!file.eof())
		{
			file >> x;
			x_s.push_back(x);
		}

		file.close();

		try
		{
			EM = new Empirical(x_s);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}

		cout << endl << endl;
		end_empirical(EM);
		break;

	case '2':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;

	}

}


void empirical_D()
{
	Empirical* ED;
	ifstream file;
	int n, k;
	char option;

	cout << "Откуда берем параметры для эмпирического распределения?" << endl;
	cout << "1. Основное распределение" << endl;
	cout << "2. Смесь распределений" << endl;
	cout << "3. На базе существующего эмпирического распределения" << endl;
	cout << "4. На базе существующей выборки элементов случайной величины" << endl;
	cout << "5. Собственные параметры" << endl;
	cout << "6. Ввести параметры из файла" << endl;
	cout << "7. Назад" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		empirical_huber();
		break;

	case '2':
		empirical_mixture();
		break;

	case '3':
		empirical_empirical();
		break;

	case '4':
		empirical_selection();
		break;

	case '5':
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		ED = new Empirical(n, k);
		end_empirical(ED);
		break;

	case '6':
		try {
			ED = new Empirical(file);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		end_empirical(ED);
		break;

	case '7':
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}