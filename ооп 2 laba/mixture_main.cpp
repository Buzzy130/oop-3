#include "empirical_main.h"

void end_mixture(Mixture* MD)
{
	ofstream x_selection;
	ofstream y_selection;
	ofstream params;
	vector<pair<double, double>> table;
	int n;
	double x = 0;

	cout << "1. ������� ��������� ����� ������������� �� �����" << endl;
	cout << "P = " << MD->get_p() << endl;
	cout << "V1 = " << MD->get_component1()->get_v() << endl;
	cout << "K1 = " << MD->get_component1()->get_k() << endl;
	cout << "Scale1 = " << MD->get_component1()->get_scale() << endl;
	cout << "Shift1 = " << MD->get_component1()->get_shift() << endl;
	cout << "V2 = " << MD->get_component2()->get_v() << endl;
	cout << "K2 = " << MD->get_component2()->get_k() << endl;
	cout << "Scale2 = " << MD->get_component2()->get_scale() << endl;
	cout << "Shift2 = " << MD->get_component2()->get_shift() << endl;
	cout << "M[X] = " << MD->Mksi_mixture() << endl;
	cout << "D[X] = " << MD->Dksi_Mixture() << endl;
	cout << "Asymmetry[X] = " << MD->asymmetry_mixture() << endl;
	cout << "Kurtosis[X] = " << MD->kurtosis_mixture() << endl;
	cout << "Density(0) = " << MD->H_Mixture(x) << endl << endl << endl;

	cout << "2. ��������� �������� ��������� ����� ������������� � ������������ �����" << endl;
	cout << "������� x: ";
	cin >> x;
	cout << endl << "Density(" << x << ") = " << MD->H_Mixture(x) << endl << endl;

	cout << "3. �������� ������� ��� �������" << endl;
	x_selection.open("selection1.txt");
	y_selection.open("selection2.txt");
	cout << "������� ����������� ������� n: ";
	cin >> n;
	table = MD->generate_pair_mixture(n);
	for (const pair<double, double>& pr : table)
	{
		x_selection << pr.first << endl;
		y_selection << pr.second << endl;
	}

	x_selection.close();
	x_selection.close();

	cout << endl << "�������� ������� �������� � ����� selection1.txt � selection2.txt" << endl << endl;
	

	cout << "4. ������� ��������� ����� ������������� � ����" << endl;
	MD->save_to_file(params);
	cout << endl << endl;


}

void mixture_D()
{
	char option;
	double p, v1, scale1, shift1, v2, scale2, shift2;
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
		end_mixture(MD);
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
		end_mixture(MD);
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
		end_mixture(MD);
		break;

	case '4':
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}