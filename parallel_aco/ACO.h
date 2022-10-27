#pragma once
#include <iostream>
#include <random>
#include <stdio.h>
#include <fstream>
#include "Livox.h"
using namespace std;
class Ant {
public:
	Ant();
	~Ant() {}
	float fitness;
	float rank;//22個變數都排50名
	vector<float>position;//22個變數
	float w;//費洛蒙濃度
	float probs;//被選中的機率	
};
class ACO {
public:
	ACO();
	~ACO() {}
	void init();
	void optimize();
	void calculateSD();
	void move();
	void compute_f();
	void update();
	void set_p_max(float* p);
	void set_p_min(float* p);
	float random(float lower, float upper);
	vector<Ant>pool;
private:
	float p_max[6];
	float p_min[6];
	int nAnts;
	int nNewAnts;
	int iteration;
	float q;//區域性權重
	float Ksi;//收斂速度因子

	float total_w;
	unsigned variable;//5+6
	vector<float>SD;
	vector<float>AVG;
	vector<float>best_position;
	vector<float>sort_fitness;
	float best_fitness;
	vector<Ant>ant;
	vector<Ant>new_ant;	
};
bool compare(int& a, int& b);
bool compare2(const Ant& a, const Ant& b);