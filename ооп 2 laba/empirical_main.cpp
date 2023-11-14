#include "empirical_main.h"


void end_empirical(Empirical* ED)
{
	ofstream file;
	double x = 0;

	cout << "������� ��������� ������������� ������������� �� �����" << endl;
	cout << "n = " << ED->get_n() << endl;
	cout << "K = " << ED->get_k() << endl;
	cout << "M[X] = " << ED->Mn() << endl;
	cout << "D[X] = " << ED->Dn() << endl;
	cout << "As[X] = " << ED->asymmetry_empirical() << endl;
	cout << "Kurt[X] = " << ED->kurtosis_empirical() << endl;
	cout << "f(0) = " << ED->H_Empirical(x) << endl << endl << endl;

	cout << "��������� �������� ��������� ����� ������������� � ������������ �����" << endl;
	cout << "������� x: ";
	cin >> x;
	cout << endl << "f(" << x << ") = " << ED->H_Empirical(x) << endl << endl;

	cout << "3. ������� ������� ��������� � ��������� ������������� ������������� � ����" << endl;
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

	cout << "�������� �����: " << endl;
	cout << "1. ����������� ������������� �������" << endl;
	cout << "2. ������ ������������ ��������� � ����������" << endl;
	cout << "3. ������� ��������� �� �����" << endl;
	cout << "4. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		cout << "������� �������� v: ";
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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;
		cout << endl << endl;

		ED = new Empirical(HD, n, k);
		end_empirical(ED);

		break;

	case '2':
		cout << "������� v scale shift: ";
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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
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

		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;

		ED = new Empirical(HD, n, k);
		end_empirical(ED);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "�������� �����: " << endl;
	cout << "1. ����������� ��������� ����� �������������" << endl;
	cout << "2. ���� � ����������" << endl;
	cout << "3. ���� �� �����" << endl;
	cout << "4. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		MD = new Mixture();

		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MD, n, k);

		end_empirical(EM);
		break;

	case '2':
		cout << "������� p (v1 scale1 shift1) (v2 scale2 shift2) (����� ������ ��� ����� Enter): ";
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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
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

		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MD, n, k);
		end_empirical(EM);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "1. ������ ��������� ������������� �� �����" << endl;
	cout << "2. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		cout << "������� ��� ����� � ����������� �������������: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("������: �� ������� ������� ����");

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
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "1. ������ ������� ��������� �� �����" << endl;
	cout << "2. �����" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "������� ��� ����� � ����������� �������������: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("������: �� ������� ������� ����");

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
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;

	}

}


void empirical_D()
{
	Empirical* ED;
	ifstream file;
	int n, k;
	char option;

	cout << "������ ����� ��������� ��� ������������� �������������?" << endl;
	cout << "1. �������� �������������" << endl;
	cout << "2. ����� �������������" << endl;
	cout << "3. �� ���� ������������� ������������� �������������" << endl;
	cout << "4. �� ���� ������������ ������� ��������� ��������� ��������" << endl;
	cout << "5. ����������� ���������" << endl;
	cout << "6. ������ ��������� �� �����" << endl;
	cout << "7. �����" << endl << endl;

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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
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
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}