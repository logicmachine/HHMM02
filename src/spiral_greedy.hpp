#pragma once
#include <vector>
#include <algorithm>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"

std::vector<int> spiral_greedy(const SparseGraph& g, const KingsGraph& kg){
	static const int dy[] = { 0, 1, 0, -1 };
	static const int dx[] = { 1, 0, -1, 0 };

	const int n = g.size(), w = kg.width(), h = kg.height();
	std::vector<int> picked(h * w), order;
	order.reserve(h * w);
	for(int i = 0, y = 0, x = 0, d = 0; i < h * w; ++i){
		picked[y * w + x] = true;
		order.push_back(y * w + x);
		const int yy = y + dy[d], xx = x + dx[d];
		if(yy < 0 || h <= yy || xx < 0 || w <= xx || picked[yy * w + xx]){
			d = (d + 1) % 4;
		}
		y += dy[d];
		x += dx[d];
	}
	std::reverse(order.begin(), order.end());

	std::vector<bool> matrix(n * n);
	for(int u = 0; u < n; ++u){
		for(const int v : g[u]){ matrix[u * n + v] = matrix[v * n + u] = true; }
	}

	std::vector<int> mapping(w * h, -1);
	std::vector<bool> used(n, false);
	mapping[order[0]] = 0;
	used[0] = true;
	for(int i = 1; i < n; ++i){
		int best_u = 0, best_delta = -1;
		for(int u = 0; u < n; ++u){
			if(used[u]){ continue; }
			int delta = 0;
			for(const int p : kg[order[i]]){
				if(mapping[p] < 0){ continue; }
				const int v = mapping[p];
				if(matrix[u * n + v]){ ++delta; }
			}
			if(delta > best_delta){
				best_delta = delta;
				best_u = u;
			}
		}
		mapping[order[i]] = best_u;
		used[best_u] = true;
	}
	return mapping;
}
