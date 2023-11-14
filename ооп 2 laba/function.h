#pragma once
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


class HuberD
{
private:

	double v;

	double k;

	double scale;

	double shift;

public:

	double get_v() const;

	void set_v(const double v);

	double get_k() const;

	double get_scale() const;

	void set_scale(const double scale);

	double get_shift() const;

	void set_shift(const double shift);

	HuberD(double v = 1, double scale = 1, double shift = 0);

	double Huber(double x) const;

	double phi(double x) const;

	double phi_lower(double x) const;

	double Mksi_huber() const;

	double Dksi_huber() const;

	double asymmetry_huber() const;

	double kurtosis_huber() const;

	double P() const;

	double K(const double v) const;

	double algorithm() const;

	HuberD(ifstream& file);

	void save_file(ofstream& file) const;

	void load_file(ifstream& file);

	vector<double> selection(const int n) const;

	vector<pair<double, double>> generate_pair(const int n, const vector<double>& x_selection = {}) const;
};









class Mixture
{
private:

	double p;

	HuberD* HD1;

	HuberD* HD2;

public:

	Mixture();

	Mixture(HuberD* HD1, HuberD* HD2, double p);

	Mixture(ifstream& file);

	HuberD* get_component1();

	HuberD* get_component2();

	double get_p() const;

	void set_p(const double p);

	void save_to_file(ofstream& file);

	double H_Mixture(const double x) const;

	double Mksi_mixture() const;

	double Dksi_Mixture() const;

	double asymmetry_mixture() const;

	double kurtosis_mixture() const;

	double algorithm_mixture() const;

	vector<double> selection_mixture(const int n) const;

	vector<pair<double, double>> generate_pair_mixture(const int n, const vector<double>& x_selection = {}) const;
};





class Empirical
{
private:

	vector<double> x_selection;

	vector<double> f_selection;

	int n = 0;

	int k = 0;

public:

	Empirical(const HuberD* HD, int _n, int _k);

	Empirical(const Mixture* MD, int _n, int _k);

	Empirical(const Empirical* ED);

	Empirical(const int _n, const int _k);

	Empirical(const vector<double>& x_selection);

	Empirical(ifstream& file);

	~Empirical();

	Empirical& operator=(const Empirical& ED);

	double random_var() const;

	vector<double> generate_x_selection() const;

	vector<double> generate_f_selection() const;

	vector<pair<double, double>> generate_pair() const;

	vector<double> get_x_selection() const;

	vector<double> get_f_selection() const;

	int get_n() const;

	int get_k() const;

	void save_to_file(ofstream& file);

	double H_Empirical(const double x_selection) const;

	double Mn() const;

	double Dn() const;

	double asymmetry_empirical() const;

	double kurtosis_empirical() const;
};




/*struct HuberDistribution
{
	double v;
	double K;
	double scale;
	double shift;
};
HuberDistribution* init_huber_distribution(double v, double K, double scale = 1, double shift = 0);
double Huber(double x, HuberDistribution* HD);
double phi(double x);
double phi_lower(double x);
double Mksi_huber(HuberDistribution* HD);
double Dksi_huber(HuberDistribution* HD);
double asymmetry_huber(HuberDistribution* HD);
double kurtosis_huber(HuberDistribution* HD);
double P(HuberDistribution* HD);
double K(double v);
double algorithm(HuberDistribution* HD);








struct Mixture
{
	double p;
	HuberDistribution* HD1;
	HuberDistribution* HD2;
};
Mixture* init_mixture(double p, double v1, double v2, double scale1 = 1, double scale2 = 1, double shift1 = 0, double shift2 = 0);
double mixture_ditribution(double x, Mixture* MD);
double 𝑀ksi_mixture(Mixture* MD);
double Dksi_mixture(Mixture* MD);
double asymmetry_mixture(Mixture* MD);
double kurtosis_mixture(Mixture* MD);



double Mn(int n, vector<double> x_s);
double Dn(int n, vector<double> x_s);
double asymmetry_empirical(int n, vector<double> x_selection);
double kurtosis_empirical(int n, vector<double> x_selection);
double empirical(int n, double x, vector<double> x_selection);
vector<double> generate_sequence(int n, HuberDistribution* HD);*/