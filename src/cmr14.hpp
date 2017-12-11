// A practical heuristic for finding graph minors
// Jun Cai, Bill Macready, Aidan Roy (2014)
#pragma once
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <limits>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"
#include "xorshift128.hpp"

namespace cmr14_detail {

static std::vector<double> make_weight_table(
	const std::vector<std::vector<int>>& vertex_models,
	const KingsGraph& kg)
{
	const int h = kg.height(), w = kg.width(), d = std::max(h, w);
	std::vector<double> weights(h * w, 1.0);
	for(const auto& vm : vertex_models){
		for(const int v : vm){ weights[v] *= d; }
	}
	return weights;
}

static std::vector<int> find_minimal_vertex_model(
	int u,
	const std::vector<std::vector<int>>& vertex_models,
	const SparseGraph& g,
	const KingsGraph& kg)
{
	using pair_type = std::pair<double, int>;
	const int n = g.size(), h = kg.height(), w = kg.width();

	bool all_empty = true;
	for(const int v : g[u]){
		const auto& vm = vertex_models[v];
		if(!vm.empty()){ all_empty = false; }
	}
	if(all_empty){
		return std::vector<int>(1, modulus_random(h * w));
	}

	const auto w_table = make_weight_table(vertex_models, kg);
	std::vector<std::vector<int>> prev_tables;
	std::vector<double> w_sum(h * w);
	for(const int v : g[u]){
		const auto& vm = vertex_models[v];
		if(vm.empty()){ continue; }

		std::vector<double> wt(h * w, std::numeric_limits<double>::infinity());
		std::priority_queue<pair_type, std::vector<pair_type>, std::greater<pair_type>> pq;
		for(const int q : vm){
			wt[q] = w_table[q];
			pq.emplace(wt[q], q);
		}

		std::vector<int> prev(h * w, -1);
		while(!pq.empty()){
			const auto cur = pq.top();
			pq.pop();
			const double c = cur.first;
			const int q = cur.second;
			if(c > wt[q]){ continue; }
			for(const int p : kg[q]){
				const double nc = c + w_table[p];
				if(nc < wt[p]){
					wt[p] = nc;
					prev[p] = q;
					pq.emplace(nc, p);
				}
			}
		}

		for(int i = 0; i < h * w; ++i){ w_sum[i] += wt[i]; }
		prev_tables.emplace_back(std::move(prev));
	}

	const int gs = std::min_element(w_sum.begin(), w_sum.end()) - w_sum.begin();
	std::vector<int> vm;
	vm.push_back(gs);
	for(const auto& prev : prev_tables){
		int cur = gs;
		while(prev[cur] >= 0){
			vm.push_back(cur);
			cur = prev[cur];
		}
	}
	std::sort(vm.begin(), vm.end());
	vm.erase(std::unique(vm.begin(), vm.end()), vm.end());
	return vm;
}

static int count_max_conflicts(
	const std::vector<std::vector<int>>& vertex_models,
	const KingsGraph& kg)
{
	const int h = kg.height(), w = kg.width(), d = std::max(h, w);
	std::vector<int> counters(h * w);
	for(const auto& vm : vertex_models){
		for(const int v : vm){ ++counters[v]; }
	}
	return *std::max_element(counters.begin(), counters.end());
}

static std::vector<int> find_minor_embedding(const SparseGraph& g, const KingsGraph& kg){
	const int n = g.size(), h = kg.height(), w = kg.width();
	const int diameter = std::max(h, w);

	// randomize the vertex order x1, ..., xn
	std::vector<int> order(n);
	std::iota(order.begin(), order.end(), 0);
	for(int i = n - 1; i > 0; --i){
		std::swap(order[i], order[modulus_random(i + 1)]);
	}

	// for i \in {1, ..., n} do
	//   set \phi(x_i) := {}
	std::vector<std::vector<int>> vertex_models(n);

	std::pair<int, int> last_score(
		std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
	bool improving = true;
	while(improving){
		for(const int i : order){
			vertex_models[i].clear();
			vertex_models[i] = find_minimal_vertex_model(i, vertex_models, g, kg);
		}
		int sum_model_sizes = 0;
		for(const auto& vm : vertex_models){ sum_model_sizes += vm.size(); }
		std::pair<int, int> new_score(
			count_max_conflicts(vertex_models, kg), sum_model_sizes);

		improving = (new_score < last_score);
		last_score = new_score;
	}
	const int final_conflicts = count_max_conflicts(vertex_models, kg);
	if(final_conflicts <= 1){
		std::vector<int> mapping(h * w, -1);
		for(int i = 0; i < n; ++i){
			for(const int q : vertex_models[i]){ mapping[q] = i; }
		}
		return mapping;
	}else{
		return {};
	}
}

}

std::vector<int> cmr14(const SparseGraph& g, const KingsGraph& kg){
	return cmr14_detail::find_minor_embedding(g, kg);
}
