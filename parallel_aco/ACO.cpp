#include "ACO.h"
#include <cmath>
#include <time.h>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <functional>
#include <omp.h>
random_device rd;
mt19937_64 generator(rd());
using namespace std;

Ant::Ant() {
	fitness = 0;
	rank = 0;//22個變數都排50名
	w = 0;
	probs = 0;
}

void ACO::set_p_max(float* p)
{
	for (unsigned i = 0; i < 6U; i++)
	{
		p_max[i] = *(p + i);
	}
}
void ACO::set_p_min(float* p)
{
	for (unsigned i = 0; i < 6U; i++)
	{
		p_min[i] = *(p + i);
	}
}
ACO::ACO() {
	nAnts =0;
	nNewAnts = 0;
	iteration = 0;
	q = 0;
	Ksi = 0;
	total_w = 0;
	variable = 0;
	best_fitness = 0;
}
void ACO::init()
{
	nAnts = 25;
	nNewAnts = 25;
	iteration = 100 ;
	q = 0.5;//0.1
	Ksi = 1;
	total_w = 0;
	variable = nMid70 * 5 + nHorizon * 6;
	best_fitness = -1;
	sort_fitness.resize(nAnts);
	SD.resize(variable);
	AVG.resize(variable);
	best_position.resize(variable);
	ant.resize(nAnts);
	new_ant.resize(nNewAnts);
	for (int i = 0; i < nAnts; i++) {
		ant[i].position.resize(variable);
	}
	for (int i = 0; i < nNewAnts; i++) {
	 	new_ant[i].position.resize(variable);
	}
	for (unsigned i = 0; i < nAnts; i++)
	{
		for (unsigned j = 0; j < 10; j++)
		{
			ant[i].position[j] = random(p_min[j % 5], p_max[j % 5]);
			ant[i].rank = 100;
		}
		for (unsigned j = 10; j < 22; j++)
		{
			ant[i].position[j] = random(p_min[(j + 2) % 6], p_max[(j + 2) % 6]);
		}
	}
	ant[0].position[0] = -3.52087;
	ant[0].position[1] = 1.41879;
	ant[0].position[2] = 0;
	ant[0].position[3] = 26.5538;
	ant[0].position[4] = -11.8823;

	ant[0].position[5] = -3.46;
	ant[0].position[6] = 1.09589;
	ant[0].position[7] = 0;
	ant[0].position[8] = -0.65602;
	ant[0].position[9] = -29.0591;

	ant[0].position[10] = -3.4482;
	ant[0].position[11] = 0.55;
	ant[0].position[12] = 0;
	ant[0].position[13] = 2.6818;
	ant[0].position[14] = -10.4729;
	ant[0].position[15] = -31.2571;

	ant[0].position[16] = -3.4688;
	ant[0].position[17] = 1.43768;
	ant[0].position[18] = 0;
	ant[0].position[19] = 1.00911;
	ant[0].position[20] = -13.6424;
	ant[0].position[21] = -9.17433;
	compute_f();
	update();
}

void ACO::optimize() {
	t = 0.0f, t_s = 0.0f, t_p = 0.0f;
	auto start_iter = chrono::steady_clock::now();
	std::fstream opt2("best_position.csv", std::ios::app | std::ios::out);
	if (!opt2.good())std::cerr << "file failed!\n";
	opt2 << nAnts << ", " << nNewAnts << "\n";
	for (unsigned i = 0; i < iteration; i++)
	{
		calculateSD();
		move();
		compute_f();
		update();
		cout << setfill('0') << setw(3) << i << setfill(' ') << setw(5) << best_fitness << "\n";
		cout << " best_position=";
		for (size_t i = 0; i < variable; i++)
		{
			cout << best_position[i] << ", ";
		}
		cout << "\n";
		opt2 << i << ", " << best_fitness << ", ";
		for (size_t i = 0; i < variable; i++)
		{
			opt2 << best_position[i] << ", ";
		}
		opt2 << "\n";
	}
	opt2.close();
	auto time_iter = chrono::duration<float>(chrono::steady_clock::now() - start_iter);
	t = time_iter.count() / iteration;
	t_p /= iteration;
	t_s = t - t_p;
	cerr << "\nt " << t << ", tp " << t_p << ", ts " << t_s << endl;

	cloud.clear();
	std::vector<Livox> Sensors;
	Sensors.reserve(nMid70 + nHorizon);
	for (float mid = 0; mid < nMid70; mid++)
	{
		Mid70 tmp_mid70(best_position.begin() + 5U * mid);
		Sensors.emplace_back(tmp_mid70);
	}
	for (unsigned horizon = 0; horizon < nHorizon; horizon++)
	{
		Horizon best02(best_position.begin() + 5U * nMid70 + 6U * horizon);
		Sensors.emplace_back(best02);
	}

	std::fstream opt("diff_amount.csv", std::ios::app | std::ios::out);
	if (!opt.good())std::cerr << "file failed!\n";
	opt  << nMid70 << ", " << nHorizon << ", \n";
	opt << best_fitness << ", ";
	for (size_t i = 0; i < variable; i++)
	{
		opt << best_position[i] << ", ";
	}
	opt << "\n";
	opt.close();

	switch (flag_cross_section)
	{
	case choose_arch:
		Arch::statistics(true);
		break;
	case choose_cylinder:
		Cylinder::statistics(true);
		break;
	case choose_horseshoe:
		Horseshoe::statistics(true);
		break;
	case choose_square:
		Square::statistics(true);
		break;
	default:
		break;
	}
	//Livox::visualiation(Sensors);
}

void ACO::compute_f()
{	
#pragma omp parallel for num_threads(THREADS) schedule(dynamic, nAnts/THREADS)
	for (unsigned i = 0; i < nAnts; i++)
	{		
		cloud.clear();
		for (unsigned mid = 0; mid < nMid70; mid++)
		{
			Mid70 tmp_mid70(ant[i].position.begin() + mid * 5U);			
			//point[mid] = nIntersections;
		}
		for (unsigned hori = 0; hori < nHorizon; hori++)
		{
			Horizon tmp_horizon(ant[i].position.begin() + nMid70 * 5U + hori * 6U);
		}
		switch (flag_cross_section)
		{
		case choose_arch:
			ant[i].fitness = Arch::statistics();
			sort_fitness[i] = ant[i].fitness;
			break;
		case choose_cylinder:
			ant[i].fitness = Cylinder::statistics();
			sort_fitness[i] = ant[i].fitness;
			break;
		case choose_horseshoe:
			ant[i].fitness = Horseshoe::statistics();
			sort_fitness[i] = ant[i].fitness;
			break;
		case choose_square:
			ant[i].fitness = Square::statistics();
			sort_fitness[i] = ant[i].fitness;
			break;
		default:
			break;
		}
	}
	sort(sort_fitness.begin(), sort_fitness.end(), greater<int>());//sort_fitness裡面放最好到最不好	

	for (int i = 0; i < nAnts; i++) {
		for (int j = 0; j < nAnts; j++) {
			if (sort_fitness[i] == ant[j].fitness) {
				ant[j].rank = i + 1;
			}
		}
	}
}

void ACO::update() {//更新費洛蒙與機率與最佳解	
	for (unsigned i = 0; i < nAnts; i++)
	{
		ant[i].w = 1 / (q * nAnts * sqrt(2 * EIGEN_PI) * exp(-(1 - ant[i].rank) * (1 - ant[i].rank) / (2 * q * q * nAnts * nAnts)));
		total_w += ant[i].w;
	}
	for (unsigned i = 0; i < nAnts; i++)
	{
		ant[i].probs = ant[i].w / total_w;
		if (ant[i].fitness > best_fitness)//這邊可以改
		{
			best_fitness = ant[i].fitness;
			for (int j = 0; j < variable; j++) {
				best_position[j] = ant[i].position[j];
			}
		}
	}
}
void ACO::calculateSD() {
	AVG.clear();
	AVG.resize(variable);
	for (unsigned i = 0; i < nAnts; i++) {
		for (unsigned j = 0; j < variable; j++)
		{
			AVG[j] += ant[i].position[j] / nAnts;
		}
	}
}

void ACO::move()
{
	double prob_sum = 0.0;
	int pos;
	// 計算機率總合
	for (int i = 0; i < nAnts; i++) {
		prob_sum += ant[i].probs;
	}
	// 計算累計機率累計分配
	for (int i = 1; i < nAnts; ++i) {
		ant[i].probs += ant[i - 1].probs;
	}
	// 開始隨機抽ant組成new_ant
	for (int i = 0; i < nNewAnts; ++i) {
		for (int j = 0; j < variable; ++j) {
			// 產生亂數
			double prob = (double)rand() / (double)RAND_MAX;
			// 取得落於哪一區塊
			for (pos = 0; pos < nAnts - 1; ++pos) {
				if (prob <= ant[pos].probs)
					break;
			}
			for (auto k = 0; k < nAnts; k++) {
				SD[j] += Ksi * abs(ant[pos].position[j] - ant[k].position[j]) / (nAnts - 1);
			}
			//新螞蟻
			if (ant[pos].position[j] == 0 && SD[j] == 0) { new_ant[i].position[j] = 0; }
			else {
				normal_distribution<float> distribution(ant[pos].position[j], SD[j]);
				new_ant[i].position[j] = (float)distribution(generator);			
				if (j < 10) {
					if (new_ant[i].position[j] > p_max[j % 5] || new_ant[i].position[j] < p_min[j % 5])
						new_ant[i].position[j] = random(p_min[j % 5], p_max[j % 5]);
				}
				else {
					if (new_ant[i].position[j] > p_max[(j + 2) % 6] || new_ant[i].position[j] < p_min[(j + 2) % 6])
						new_ant[i].position[j] = random(p_min[(j + 2) % 6], p_max[(j + 2) % 6]);
				}			
				SD[j] = 0;
			}
		}
	}
#pragma omp parallel for num_threads(THREADS) schedule(dynamic, nNewAnts/THREADS)
	for (unsigned i = 0; i < nNewAnts; i++)
	{
		cloud.clear();
		for (unsigned mid = 0; mid < nMid70; mid++)
		{
			Mid70 tmp(new_ant[i].position.begin() + mid * 5U);
		}
		for (unsigned hori = 0; hori < nHorizon; hori++)
		{
			Horizon tmp(new_ant[i].position.begin() + nMid70 * 5U + hori * 6U);
		}
		switch (flag_cross_section)
		{
		case choose_arch:
			new_ant[i].fitness = Arch::statistics();
			break;
		case choose_cylinder:
			new_ant[i].fitness = Cylinder::statistics();
			break;
		case choose_horseshoe:
			new_ant[i].fitness = Horseshoe::statistics();
			break;
		case choose_square:
			new_ant[i].fitness = Square::statistics();
			break;
		default:
			break;
		}
	}
	for (int i = 0; i < nNewAnts; ++i) {
		pool.push_back(new_ant[i]);
	}
	for (int i = 0; i < nAnts; i++) {
		pool.push_back(ant[i]);
	}	
	sort(pool.begin(), pool.end(), compare2);//全部排序
	for (int i = 0; i < nAnts; i++) {//前50名塞回去
		ant[i] = pool[i];
	}
	pool.clear();
}

float ACO::random(float lower, float upper)
{
	uniform_real_distribution<double> unif(lower, upper);
	return floor(unif(generator) * 100.0) / 100.0;
}

bool compare(int& a, int& b)
{
	return a > b;
}

bool compare2(const Ant& a, const Ant& b)
{
	return a.fitness > b.fitness;
}