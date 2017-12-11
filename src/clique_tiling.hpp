#pragma once
#include <array>
#include <vector>
#include <algorithm>
#include <chrono>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"
#include "xorshift128.hpp"

static int compute_clique_tile_size(const SparseGraph& g, const KingsGraph& kg){
	const int n = g.size(), h = kg.height(), w = kg.width();
	for(int s = n; s > 0; --s){
		const int ch = h / s, cw = h / s;
		if(ch * cw * s >= n){ return s; }
	}
	return 0;
}

static std::vector<int> clique_tiling(const SparseGraph& g, const KingsGraph& kg){
	const int n = g.size(), h = kg.height(), w = kg.width();
	const int ts = compute_clique_tile_size(g, kg), th = h / ts, tw = w / ts;

	std::vector<bool> matrix(n * n);
	for(int u = 0; u < n; ++u){
		for(const int v : g[u]){
			matrix[u * n + v] = matrix[v * n + u] = true;
		}
	}

	const int m = ts * th * tw;
	std::vector<int> emb_mapping(w * h, -1);
	for(int ty = 0; ty < th; ++ty){
		for(int tx = 0; tx < tw; ++tx){
			const int y_offset = ty * ts, x_offset = tx * ts;
			const int tile_id = ty * tw + tx;
			for(int i = 0; i < ts; ++i){
				const int u = tile_id * ts + i;
				int c = i, d = (i % 2 == 0 ? 1 : -1);
				for(int j = 0; j < ts; ++j){
					const int q = ((ty + tx) % 2 == 0)
						? ((y_offset + c) * w + (x_offset + j))
						: ((y_offset + j) * w + (x_offset + c));
					emb_mapping[q] = u;
					if(c == 0 && d < 0){
						d = -d;
					}else if(c == ts - 1 && d > 0){
						d = -d;
					}else{
						c += d;
					}
				}
			}
		}
	}

	std::vector<std::vector<int>> emb_adjacency(m);
	for(int q = 0; q < h * w; ++q){
		const int u = emb_mapping[q];
		if(u < 0){ continue; }
		for(const int p : kg[q]){
			const int v = emb_mapping[p];
			if(v < 0 || v == u){ continue; }
			emb_adjacency[u].push_back(v);
		}
	}
	for(auto& v : emb_adjacency){
		std::sort(v.begin(), v.end());
		v.erase(std::unique(v.begin(), v.end()), v.end());
	}

	std::default_random_engine engine;
	std::uniform_int_distribution<int> n_dist(0, n - 1);

	std::vector<int> mapping(m, -1);
	std::vector<int> imapping(n);
	std::vector<int> init_order(m);
	std::iota(init_order.begin(), init_order.end(), 0);
	std::shuffle(init_order.begin(), init_order.end(), engine);
	for(int i = 0; i < n; ++i){
		mapping[init_order[i]] = i;
		imapping[i] = init_order[i];
	}

	int score = 0;
	for(int q = 0; q < m; ++q){
		const int u = mapping[q];
		if(u < 0){ continue; }
		for(const int p : emb_adjacency[q]){
			const int v = mapping[p];
			if(v < 0){ continue; }
			if(matrix[u * n + v]){ ++score; }
		}
	}
	score /= 2;

	std::vector<int> best_mapping = mapping;
	int best_score = score;

	const auto time_limit = std::chrono::milliseconds(3000);
	const auto start_time = std::chrono::steady_clock::now();
	auto last_time = start_time;
	auto max_duration = last_time - start_time;

	const double max_temp = 5.0;
	const double min_temp = 2.0;

	do {
		const double progress = (last_time - start_time) / time_limit;
		const double temperature = max_temp - (max_temp - min_temp) * progress;

		constexpr int diff_limit = 64;
		std::array<double, diff_limit + 1> exp_table;
		const double e = exp(-1.0 / temperature);
		double exp_acc = 1.0;
		for(int i = 0; i <= diff_limit; ++i){
			exp_table[i] = exp_acc;
			exp_acc *= e;
		}

		for(int iter = 0; iter < 1000; ++iter){
			const int u = modulus_random(n);
			if(g[u].empty()){ continue; }
			const int v = g[u][modulus_random(g[u].size())];
			const auto& q_neighbors = emb_adjacency[imapping[u]];
			const int d = modulus_random(q_neighbors.size());
			const int q = q_neighbors[d], p = imapping[v];
			const int w = mapping[q];
			if(p == q){ continue; }
			int diff = 0;
			if(w >= 0){ // remove w from q
				for(const int r : emb_adjacency[q]){
					if(mapping[r] < 0){ continue; }
					if(matrix[w * n + mapping[r]]){ --diff; }
				}
				mapping[q] = imapping[w] = -1;
			}
			if(v >= 0){ // remove v from p
				for(const int r : emb_adjacency[p]){
					if(mapping[r] < 0){ continue; }
					if(matrix[v * n + mapping[r]]){ --diff; }
				}
				mapping[p] = imapping[v] = -1;
			}
			if(v >= 0){ // put v to q
				for(const int r : emb_adjacency[q]){
					if(mapping[r] < 0){ continue; }
					if(matrix[v * n + mapping[r]]){ ++diff; }
				}
				mapping[q] = v;
				imapping[v] = q;
			}
			if(w >= 0){ // put w to p
				for(const int r : emb_adjacency[p]){
					if(mapping[r] < 0){ continue; }
					if(matrix[w * n + mapping[r]]){ ++diff; }
				}
				mapping[p] = w;
				imapping[w] = p;
			}
			const double r = xorshift128() / static_cast<double>(0xffffffffu);
			bool undo =
				(diff < -diff_limit) || (r > exp_table[std::max(-diff, 0)]);
			if(undo){
				mapping[q] = w;
				mapping[p] = v;
				if(w >= 0){ imapping[w] = q; }
				if(v >= 0){ imapping[v] = p; }
			}else{
				score += diff;
				if(score > best_score){
					std::cerr << best_score << " -> " << score << " (" << score - best_score << ")" << std::endl;
					best_score = score;
					best_mapping = mapping;
				}
			}
		}

		const auto now = std::chrono::steady_clock::now();
		max_duration = std::max(max_duration, now - last_time);
		last_time = now;
	} while((last_time - start_time) + max_duration * 1.0 < time_limit);

	for(auto& x : emb_mapping){
		if(x >= 0){ x = mapping[x]; }
	}
	return emb_mapping;
}
