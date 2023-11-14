#include "function.h"

//--------------------------------------------------------------------------huber_distribution--------------------------------------------------------------------------------//
HuberD::HuberD(double v, double scale, double shift)
{
	if (v <= 0 || scale <= 0)
		throw invalid_argument("Ошибка: один или несколько параметров не корректны");

	this->v = v;
	this->k = K(v);
	this->scale = scale;
	this->shift = shift;
}

HuberD::HuberD(ifstream& file)
{
	load_file(file);
}

double HuberD::get_v() const
{
	return this->v;
}

double HuberD::get_k() const
{
	return this->k;
}

double HuberD::get_scale() const
{
	return this->scale;
}

double HuberD::get_shift() const
{
	return this->shift;
}

void HuberD::set_v(const double v)
{
	if (v <= 0)
		throw invalid_argument("Ошибка: v <= 0");

	this->v = v;
	this->k = K(this->v);
}

void HuberD::set_scale(const double scale)
{
	if (scale <= 0)
		throw invalid_argument("Ошибка: λ <= 0");
	this->scale = scale;
}

void HuberD::set_shift(const double shift)
{
	this->shift = shift;
}

double HuberD::Huber(const double x) const
{
	if (abs((x - this->shift) / this->scale) <= this->v)
	{
		return (1. / (sqrt(2. * M_PI) * this->k) * exp(-pow((x - this->shift) / this->scale, 2.) / 2.)) / this->scale;
	}
	if (abs((x - this->shift) / this->scale) > this->v)
	{
		return (1. / (sqrt(2. * M_PI) * this->k) * exp(pow(this->v, 2.) / 2. - this->v * abs((x - this->shift) / this->scale))) / this->scale;
	}
}
double HuberD::phi(double x) const//Ф(x)
{
	return 0.5 * (1. + erf(x / sqrt(2.)));
}
double HuberD::phi_lower(double x) const//ф(x)
{
	return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}
double HuberD::Mksi_huber() const//мат ожидание
{
	return this->shift;
}
double HuberD::Dksi_huber() const//дисперсия
{
	return 1. + 2. * phi_lower(this->v) * (pow(this->v, 2.) + 2.) / (pow(this->v, 3.) * this->k);
}
double HuberD::asymmetry_huber() const//ассиметрия
{
	return 0.;
}
double HuberD::kurtosis_huber() const//коэфф эксцесса
{
	return (3. * (2. * phi(this->v) - 1.) + 2. * phi_lower(this->v) * (24. / pow(this->v, 5.) + 24. / pow(this->v, 3.) + 12. /
		this->v + this->v)) / (pow(Dksi_huber(), 2.) * this->k) - 3.;
}
double HuberD::P() const//вероятности попадания в центральный интервал 
{
	return (2. * phi(this->v) - 1.) / this->k;
}
double HuberD::K(double v) const//значение K зная V
{
	return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}
//-----------------------------------------------------//
double HuberD::algorithm() const
{
	std::random_device rd;
	std::default_random_engine gen(rd());
	std::uniform_real_distribution<> d(0, 1);
	//шаг 1
	double r1 = d(gen);
	if (r1 <= P())
	{
		//шаг 2
		double r2, r3, x1;
		do {
			r2 = d(gen);
			r3 = d(gen);
			x1 = sqrt(-2 * log(r2)) * cos(2 * M_PI * r3);
			//double x1 = sqrt(-2 * log(r2)) * sin(2 * M_PI * r3)
		} while (!(-this->v <= x1 && x1 <= this->v)); //шаг 3
		return x1 * this->scale + this->shift;
	}
	else
	{
		//шаг 4
		double r4 = d(gen);
		double x2 = this->v - log(r4) / this->v;
		//шаг 5
		return r1 < (1 + P()) / 2 ? x2 * this->scale + this->shift : -x2 * this->scale + this->shift;
	}
}

void HuberD::load_file(ifstream& file)
{
	string filename;
	double v, scale, shift;

	file.open("input.txt");

	if (!file)
		throw runtime_error("Ошибка: не удалось открыть файл");

	file >> v >> scale >> shift;

	file.close();

	if (v <= 0 || scale <= 0)
		throw invalid_argument("Ошибка: один или несколько параметров некорректны");

	this->v = v;
	this->scale = scale;
	this->shift = shift;
}

void HuberD::save_file(ofstream& file) const
{
	string filename = "output.txt";

	file.open(filename);
	file << this->v << endl << this->scale << endl << this->shift;
	file.close();
}

vector<double> HuberD::selection(const int n) const
{
	vector<double> res;

	for (int i = 0; i < n; i++)
	{
		double x = algorithm();
		res.push_back(x);
	}

	sort(res.begin(), res.end());

	return res;
}

vector<pair<double, double>> HuberD::generate_pair(const int n, const vector<double>& x_selection) const
{
	vector<pair<double, double>> res;
	vector<double> sequence;

	if (x_selection.empty())//если строка пуста заполни ее
		sequence = selection(n);
	else
		sequence = x_selection;

	for (const double& x : sequence)
	{
		double y = Huber(x);
		res.push_back(make_pair(x, y));
	}

	return res;
}





//--------------------------------------------------------------------------mixture_distribution--------------------------------------------------------------------------------//


Mixture::Mixture()
{
	p = 0.5;
	HD1 = new HuberD(1, 1, 0);
	HD2 = new HuberD(1, 1, 0);
}

Mixture::Mixture(HuberD* _HB1, HuberD* _HB2, double _p) :
	HD1(_HB1), HD2(_HB2), p(_p > 0 and _p < 1 ? _p : throw invalid_argument("Ошибка: один или несколько параметров некорректны")) {}

Mixture::Mixture(ifstream& file)
{
	double v1, scale1, shift1, v2, scale2, shift2, p;
	string filename = "input.txt";
	file.open(filename);

	if (!file)
		throw runtime_error("Ошибка: не удалось открыть файл");

	file >> p >> v1 >> scale1 >> shift1 >> v2 >> scale2 >> shift2;
	file.close();

	if (p > 1 || p < 0 || v1 <= 0 || v2 <= 0 || scale1 <= 0 || scale2 <= 0)
		throw invalid_argument("Ошибка: один или несколько параметров не корректны");

	HD1 = new HuberD(v1, scale1, shift1);
	HD2 = new HuberD(v2, scale2, shift2);

	this->p = p;
}

HuberD* Mixture::get_component1()
{
	return HD1;
}

HuberD* Mixture::get_component2()
{
	return HD2;
}

double Mixture::get_p() const
{
	return p;
}

void Mixture::set_p(const double p)
{
	if (p > 1 || p < 0)
		throw invalid_argument("Ошибка: один параметро некорректен");

	this->p = p;
}




double Mixture::H_Mixture(const double x) const//плотность
{
	return (1 - p) * HD1->Huber(x) + p * HD2->Huber(x);
}

double Mixture::Mksi_mixture() const
{
	return (1 - p) * HD1->Mksi_huber() + p * HD2->Mksi_huber();
}

double Mixture::Dksi_Mixture() const
{
	return (1 - p) * (pow(HD1->Mksi_huber(), 2) + HD1->Dksi_huber()) +
		p * (pow(HD2->Mksi_huber(), 2) + HD2->Dksi_huber()) -
		pow(Mksi_mixture(), 2);


}

double Mixture::asymmetry_mixture() const
{
	return (1 / pow(Dksi_Mixture(), 3 / 2)) *
		((1 - p) *
			(pow(HD1->Mksi_huber() - Mksi_mixture(), 3) + 3 *
				(HD1->Mksi_huber() - Mksi_mixture()) * HD1->Dksi_huber() +
				pow(HD1->Dksi_huber(), 3 / 2) * HD1->asymmetry_huber()) +
			(p) *
			(pow(HD2->Mksi_huber() - Mksi_mixture(), 3) + 3 *
				(HD2->Mksi_huber() - Mksi_mixture()) * HD2->Dksi_huber() +
				pow(HD2->Dksi_huber(), 3 / 2) * HD2->asymmetry_huber()
				));



}

double Mixture::kurtosis_mixture() const
{
	return (1 / pow(Dksi_Mixture(), 2)) * (
		(1 - p) * (pow(HD1->Mksi_huber() - Mksi_mixture(), 4) +
			6 * pow(HD1->Mksi_huber() - Mksi_mixture(), 2) * HD1->Dksi_huber() +
			4 * (HD1->Mksi_huber() - Mksi_mixture()) * pow(HD1->Dksi_huber(), 3 / 2) * HD1->asymmetry_huber() +
			pow(HD1->Dksi_huber(), 2) * HD1->kurtosis_huber()) +
		(p) * (pow(HD2->Mksi_huber() - Mksi_mixture(), 4) +
			6 * pow(HD2->Mksi_huber() - Mksi_mixture(), 2) * HD2->Dksi_huber() +
			4 * (HD2->Mksi_huber() - Mksi_mixture()) * pow(HD2->Dksi_huber(), 3 / 2) * HD2->asymmetry_huber() +
			pow(HD2->Dksi_huber(), 2) * HD2->kurtosis_huber()) - 3);
}



void Mixture::save_to_file(ofstream& file)
{
	string filename = "output.txt";
	file.open(filename);

	file << p << endl << HD1->get_v() << endl << HD1->get_scale() << endl << HD1->get_shift() << endl << HD2->get_v() << endl << HD2->get_scale() << endl << HD2->get_shift();
	file.close();

	cout << endl << "Параметры распределения сохранены в файл " + filename << endl;
}


double Mixture::algorithm_mixture() const {
	random_device rd;
	default_random_engine gen(rd());
	uniform_real_distribution<> d(0, 1);

	double r = d(gen);

	if (r > p)
		return HD1->algorithm();
	else
		return HD2->algorithm();
}

vector<double> Mixture::selection_mixture(const int n) const
{
	vector<double> sequence;
	for (int i = 0; i < n; i++)
	{
		double x = algorithm_mixture();
		sequence.push_back(x);
	}

	sort(sequence.begin(), sequence.end());
	return sequence;
}

vector<pair<double, double>> Mixture::generate_pair_mixture(const int n, const vector<double>& x_selection) const
{
	vector<double> sequence;
	vector<pair<double, double>> table;

	if (x_selection.empty())
		sequence = selection_mixture(n);
	else
		sequence = x_selection;

	for (const double& x : sequence)
	{
		double y = H_Mixture(x);
		table.push_back(make_pair(x, y));
	}

	return table;
}

//--------------------------------------------------------------------------empirical_distribution--------------------------------------------------------------------------------//

vector<double> Empirical::generate_x_selection() const//Функция генерирует и возвращает вектор result, содержащий случайные числа.
{
	vector<double> result;

	for (int i = 0; i < n; i++)
		result.push_back(random_var());//от 0 до n-1 генерируем числа типа double и добавляется в вектор result

	sort(result.begin(), result.end());//сортируем в порядке возрастания

	return result;
}


vector<double> Empirical::generate_f_selection() const//генерирование значений плотности по выборке элеметов случайной величины для H_Empirical
{
	vector<double> result;//Создается пустой вектор result, который будет содержать значения функции H_Empirical().

	for (const double& x : x_selection)
		result.push_back(H_Empirical(x));
	//После завершения цикла, вектор result содержит значения функции H_Empirical() для каждого элемента x в векторе x_selection.

	return result;
}

Empirical::Empirical(const HuberD* HD, int _n, int _k)//генерирование выборки по хьюберу
{
	if (_n <= 1) {
		throw invalid_argument("Некорректный аргумент");
	}
	if (_k <= 2) {//Если _k меньше или равно 2, то значение _k пересчитывается с использованием формулы _k = static_cast<int>(floor(log2(_n))) + 1
		_k = static_cast<int>(floor(log2(_n))) + 1;
	}

	n = _n;//размер выборки
	k = _k;//кол-во интервалов
	x_selection = HD->selection(_n);//Этот метод генерирует и возвращает вектор x_selection
	f_selection = generate_f_selection();//Этот вектор инициализирует член данных f_selection в классе Empirical
}

Empirical::Empirical(const Mixture* HD, int _n, int _k)//генерирование выборки по смеси
{
	if (_n <= 1) {
		throw invalid_argument("Некорректный аргумент");
	}

	n = _n;
	k = (_k > 2) ? _k : static_cast<int>(floor(log2(_n))) + 1;//Если _k больше 2, то значение _k присваивается члену данных k в классе Empirical. Если _k меньше или равно 2, то значение _k пересчитывается с использованием формулы _k = static_cast<int>(floor(log2(_n))) + 1
	x_selection = HD->selection_mixture(_n);//Этот метод генерирует и возвращает вектор x_selection
	f_selection = generate_f_selection();//Этот вектор инициализирует член данных f_selection в классе Empirical
}

Empirical::Empirical(const Empirical* ED)//конструктор копирования
{
	if (ED->n <= 1) {
		throw invalid_argument("Некорректный аргумент");
	}

	n = ED->n;//n из объекта ED присваивается члену данных n в классе Empirical
	k = (ED->k > 2) ? ED->k : static_cast<int>(floor(log2(ED->n))) + 1;//Если значение k больше 2, то значение k присваивается члену данных k в классе Empirical. Если значение k меньше или равно 2, то значение k пересчитывается с использованием формулы k = static_cast<int>(floor(log2(ED->n))) + 1
	x_selection = ED->x_selection;//x_selection из объекта ED присваивается члену данных x_selection в классе Empirical
	f_selection = generate_f_selection();//генерирует вектор значений функции H_Empirical() для каждого элемента вектора x_selection. Этот вектор инициализирует член данных f_selection в классе Empirical.
}

Empirical::Empirical(const int _n, const int _k)//генерирование выборки по эмпирическому распределению
{
	if (_n <= 1) {
		throw std::invalid_argument("Некорректный аргумент");
	}

	n = _n;
	k = (_k > 2) ? _k : static_cast<int>(std::floor(std::log2(_n))) + 1;//если K>2 k=_k иначе k = static_cast<int>(std::floor(std::log2(_n))) + 1
	x_selection = generate_x_selection();
	f_selection = generate_f_selection();
}

Empirical::Empirical(const vector<double>& x_selection)//конструктор по выборке
	: n(x_selection.size()), k(static_cast<int>(std::floor(std::log2(x_selection.size()))) + 1)//n=x_selection.size() k=static_cast<int>(std::floor(std::log2(x_selection.size()))
{
	this->x_selection = x_selection;
	f_selection = generate_f_selection();
}

Empirical::~Empirical()//деструктор
{
	x_selection.clear();
	f_selection.clear();
}

Empirical& Empirical::operator=(const Empirical& ED)//перегрузка оператора присваивания
{
	if (this == &ED)//является ли текующий обьект и переданный одинаковыми
		return *this;//возвращаем указатель на текующий обьект
	//иначе присваиваем текущим обьектам
	x_selection = ED.x_selection;
	f_selection = ED.f_selection;
	n = ED.n;
	k = ED.k;
	return *this;
	//Таким образом, этот оператор присваивания позволяет скопировать значения членов данных из одного объекта Empirical в другой объект Empirical того же типа.
}

double Empirical::random_var() const//генерирует случайную величину на основе эмпирического распределения.
{
	vector<double> intervals;//содержать интервалы, на которые будет разбит отрезок [0, 1) для генерации случайной величины.
	vector<double> densities;//будет содержать плотности вероятности для каждого интервала.

	int k = get_k();//возвращает количество интервалов.
	double delta = (1.0 - 0.0) / (double)k;// которое представляет ширину каждого интервала на отрезке [0, 1).
	//Затем происходит генерация интервалов
	double point = 0.0;//Изначально устанавливается значение point равное 0.0.
	//добавляются значения point в вектор intervals и увеличивается на значение delta, пока point не превышает 1.0.
	while (point < 1.0)
	{
		intervals.push_back(point);
		point += delta;
	}
	//intervals.push_back(1.0);

	for (int i = 0; i < k; i++)
		densities.push_back(delta);//осле генерации интервалов, создается вектор densities в котором каждая плотность вероятности равна delta.

	random_device rd;
	default_random_engine gen(rd());
	piecewise_constant_distribution<> d(intervals.begin(), intervals.end() - 1, densities.begin());

	return d(gen);

}


vector<pair<double, double>> Empirical::generate_pair() const//Этот метод используется для генерации таблицы значений, представленной в виде вектора пар (x, f(x)).
{
	vector<pair<double, double >> result;//тут храним пары значений

	for (int i = 0; i < n; i++)
		result.push_back(make_pair(x_selection[i], f_selection[i]));//добавляем каждую пару в вектор result

	return result;
}




vector<double> Empirical::get_x_selection() const
{
	return x_selection;
}

vector<double> Empirical::get_f_selection() const
{
	return f_selection;
}

int Empirical::get_n() const
{
	return n;
}

int Empirical::get_k() const
{
	return k;
}


double Empirical::H_Empirical(const double x) const//Этот метод используется для вычисления эмпирической функции распределения для данного значения x. типо плотность
{//Вычисление количества интервалов k
	int k = (int)trunc(log2((double)n)) + 1;//кол-во промежутков
	double x_min = *min_element(x_selection.begin(), x_selection.end());
	double x_max = *max_element(x_selection.begin(), x_selection.end());
	double delta = (1 / (double)k) * (x_max - x_min);//ширина каждого интервала
	//проверяется, в каком интервале находится значение x. Если значение x попадает в интервал (x_min + delta * i, x_min + delta * (i + 1)), то выполняются следующие действия:
	for (int i = 0; i < k; i++)
		if (x_min + delta * i <= x && x < x_min + delta * (i + 1))
		{
			int n_i = 0;//Инициализируется счетчик n_i для подсчета количества элементов в данном интервале.
			for (double x : x_selection) {//Во внутреннем цикле for, проходящем по каждому элементу x в векторе x_selection, проверяется, принадлежит ли элемент интервалу (x_min + delta * i, x_min + delta * (i + 1)). Если да, величивается счетчик n_i на 1.
				if (i == k - 1) {
					if (x_min + delta * (double)i <= x && x <= x_min + delta * (double)(i + 1)) {
						n_i++;
					}
				}
				else {
					if (x_min + delta * (double)i <= x && x < x_min + delta * (double)(i + 1)) {
						n_i++;
					}
				}
			}
			return n_i / (n * delta);//Возвращается значение n_i / (n * delta), которое представляет относительную частоту элементов в данном интервале.
		}
	x_selection;
	return 0.0;//Если значение x не попадает ни в один из интервалов, возвращается значение 0.0.
}

double Empirical::Mn() const
{
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += x_selection[i];
	}
	return sum / (double)n;
}


double Empirical::Dn() const//+
{
	const double Mn_ksi = Mn();
	double sum = 0;

	for (int i = 0; i < n; ++i) {
		double diff = x_selection[i] - Mn_ksi;
		sum += diff * diff;
	}

	return sum / n;
}

double Empirical::asymmetry_empirical() const
{
	const double Mn_ksi = Mn();
	const double Dn_ksi = Dn();
	double sum = 0;

	for (int i = 0; i < n; ++i)
		sum += pow(x_selection[i] - Mn_ksi, 3);

	return sum / (n * pow(Dn_ksi, 3 / 2));
}

double Empirical::kurtosis_empirical() const
{
	const double Mn_ksi = Mn();
	const double Dn_ksi = Dn();
	double sum = 0;

	for (int i = 0; i < n; ++i) {
		sum += pow(x_selection[i] - Mn_ksi, 4);
	}
	return (sum / (n * pow(Dn_ksi, 2))) - 3;
}

Empirical::Empirical(ifstream& file)
{
	string filename = "input.txt";

	file.open(filename);
	if (!file)
		throw runtime_error("Ошибка: не удалось открыть файл");

	double x;
	while (!file.eof())
	{
		file >> x;
		x_selection.push_back(x);
	}
	file.close();

	n = x_selection.size();
	k = (int)floor(log2(n)) + 1;
	f_selection = generate_f_selection();
}

void Empirical::save_to_file(ofstream& file)
{
	vector<pair<double, double>> pairs = generate_pair();

	ofstream file_selection1;
	ofstream file_selection2;
	ofstream params;

	file_selection1.open("selection1.txt");
	file_selection2.open("selection2.txt");
	params.open("params.txt");

	for (const pair<double, double>& pair : pairs)
	{
		file_selection1 << pair.first << endl;
		file_selection2 << pair.second << endl;
	}

	params << n << k;

	file_selection1.close();
	file_selection2.close();
	params.close();

	cout << "Выборка сохранена в файл selection1.txt" << endl;
	cout << "Значения плотности сохранены в selection2.txt" << endl;
	cout << "Параметры распределения сохранены в файл params.txt" << endl << endl;
}


