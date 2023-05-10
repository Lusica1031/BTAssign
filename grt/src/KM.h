#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

class KM {
public:
	//	KM(float* data, int m, int n);
	KM(std::vector<std::vector<int>> data, int m, int n) {
		init(data, m, n);
	}

	~KM();

	int N;
	int front;
	int back;
	int* matchX;
	int* matchY;
	std::vector<std::vector<int>> weights;

	void init(std::vector<std::vector<int>> data, int m, int n);
	void del();
	void compute();
	int maxWeight() {
		return max_w;
	}
	vector<int> getMatch(bool front2back = true);

private:
	
	int max_w;	
	int* flagX;
	int* flagY;
	char* usedX;
	char* usedY;
	
	void constructMatrix(std::vector<std::vector<int>> data, int m, int n);
	bool dfs(int v);
};