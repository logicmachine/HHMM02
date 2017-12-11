#pragma once
#include <vector>
#include <algorithm>
#include <memory>
#include <cstdint>
#include <cassert>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"
#include "matrix.hpp"
#include "radix_sort.hpp"

namespace skew_cyclic_detail {

struct Node {
	int score;
	std::bitset<512> status;
	std::vector<int> history;

	Node()
		: score(0)
		, status()
		, history()
	{ }
};

using NodePtr = std::shared_ptr<Node>;

struct NodeComparator {
	bool operator()(const NodePtr& a, const NodePtr& b){
		return a->score < b->score;
	}
};

static std::vector<int> make_next_table(const KingsGraph& kg){
	const int m = kg.height();
	std::vector<int> next_table(m * m, -1);
	for(int y = 0; y + 1 < m; ++y){
		for(int x = 0; x < m; ++x){
			const int q = y * m + x;
			if((y + x) % 2 == 1){
				next_table[q] = (x == m - 1)
					? ((y + 1) * m + x)
					: ((y + 1) * m + (x + 1));
			}else{
				next_table[q] = (x == 0)
					? ((y + 1) * m + x)
					: ((y + 1) * m + (x - 1));
			}
		}
	}
	return next_table;
}

static std::vector<int> make_segment_mapping(int n, const KingsGraph& kg){
	const int m = kg.height();
	const auto next_table = make_next_table(kg);
	std::vector<int> order;
	for(int i = (m - 1) / 2; i >= 0; --i){
		order.push_back(i);
		if(m - 1 - i != i){ order.push_back(m - 1 - i); }
	}
	std::vector<int> mapping(m * m, -1);
	for(int i = 0, v = 0; i < m; ++i){
		const int n_ranges = (n + m - 1) / m;
		int cur = order[i];
		for(int j = 0; j < n_ranges; ++j){
			const int y = m / n_ranges + (j < m % n_ranges ? 1 : 0);
			for(int c = 0; c < y; ++c){
				mapping[cur] = v;
				cur = next_table[cur];
			}
			if(n % m == 0 || i < n % m || j != (n_ranges - 1) / 2){ ++v; }
		}
	}
	return mapping;
}

static SparseGraph make_segment_adjacency_list(
	const std::vector<int>& mapping,
	const KingsGraph& kg)
{
	const int n = 1 + *std::max_element(mapping.begin(), mapping.end());
	Matrix<int> adj_matrix(n, n);
	for(int q = 0; q < kg.size(); ++q){
		const int u = mapping[q];
		if(u < 0){ continue; }
		for(const int p : kg[q]){
			const int v = mapping[p];
			if(v < 0 || u == v){ continue; }
			adj_matrix(u, v) = 1;
		}
	}
	SparseGraph adj_list(n);
	for(int u = 0; u < n; ++u){
		for(int v = 0; v < n; ++v){
			if(adj_matrix(u, v)){ adj_list.add_edge(u, v); }
		}
	}
	return adj_list;
}

static std::vector<int> compute_placement_order(
	const SparseGraph& g_emb,
	const std::vector<int>& mapping)
{
	const int n = g_emb.size(), m = static_cast<int>(sqrt(mapping.size()));
	const int first = mapping[(m / 2) * (m + 1)];
	std::vector<int> order;
	order.reserve(n);
	std::vector<int> degrees(n, 0);
	order.push_back(first);
	degrees[first] = -1;
	for(const int p : g_emb[first]){
		if(degrees[p] >= 0){ ++degrees[p]; }
	}
	for(int i = 1; i < n; ++i){
		const int q = max_element(degrees.begin(), degrees.end()) - degrees.begin();
		order.push_back(q);
		degrees[q] = -1;
		for(const int p : g_emb[q]){
			if(degrees[p] >= 0){ ++degrees[p]; }
		}
	}
	return order;
}

static std::vector<int> inverse_permutation(const std::vector<int>& p){
	const int n = p.size();
	std::vector<int> q(n);
	for(int i = 0; i < n; ++i){ q[p[i]] = i; }
	return q;
}

std::vector<int> run(const SparseGraph& g, const KingsGraph& kg){
	assert(kg.width() == kg.height());
	const int n = g.size(), m = kg.width();
	const auto segment_mapping = make_segment_mapping(n, kg);
	const auto g_emb = make_segment_adjacency_list(segment_mapping, kg);
	const auto order = compute_placement_order(g_emb, segment_mapping);
	const auto inv_order = inverse_permutation(order);

	Matrix<uint8_t> gmat(n, n);
	for(int u = 0; u < n; ++u){
		for(const int v : g[u]){ gmat(u, v) = 1; }
	}

	// beam search
	static const int beam_width = 4000;
	RadixSortEngine<int, std::pair<int, int>> radix_sort;
	std::vector<Node> states(1);
	for(int i = 0; i < n; ++i){
		const int q = order[i];
		const auto& neighbors = g_emb[q];

		std::vector<int> next_scores;
		std::vector<std::pair<int, int>> next_moves;
		next_scores.reserve(states.size() * (n - i));
		next_moves.reserve(states.size() * (n - i));

		for(size_t j = 0; j < states.size(); ++j){
			const auto& cur = states[j];
			for(int u = 0; u < n; ++u){
				if(cur.status[u]){ continue; }
				int s = cur.score;
				for(const int p : neighbors){
					if(inv_order[p] >= i){ continue; }
					const int v = cur.history[inv_order[p]];
					s += gmat(u, v);
				}
				next_scores.push_back(s);
				next_moves.emplace_back(j, u);
			}
		}
		const size_t nn = next_scores.size();
		radix_sort(nn, next_scores.data(), next_moves.data());
		const size_t lo = nn > beam_width ? nn - beam_width : 0, hi = nn;
		std::vector<Node> next_states;
		next_states.reserve(hi - lo);
		for(size_t s = lo; s < hi; ++s){
			const int score = next_scores[s];
			const int j = next_moves[s].first;
			const int u = next_moves[s].second;
			const auto& cur = states[j];
			Node next;
			next.score = score;
			next.status = cur.status;
			next.status[u] = true;
			next.history.reserve(i + 1);
			for(const auto& x : cur.history){ next.history.push_back(x); }
			next.history.push_back(u);
			next_states.push_back(std::move(next));
		}

		states = std::move(next_states);
	}

	const int score = states.back().score;
	std::vector<int> mapping(n), imapping(n);
	for(int i = 0; i < n; ++i){
		const int v = states.back().history[i];
		mapping[v] = order[i];
		imapping[order[i]] = v;
	}

	int score_limit_emb = 0, score_limit = 0;
	for(int u = 0; u < n; ++u){
		score_limit += g[u].size();
		score_limit_emb += g_emb[u].size();
	}
	score_limit /= 2;
	score_limit_emb /= 2;
	std::cerr << score << " / (" << score_limit << ", " << score_limit_emb << ")" << std::endl;

	std::vector<int> solution(m * m, -1);
	for(int q = 0; q < m * m; ++q){
		const int s = segment_mapping[q];
		if(s < 0){ continue; }
		solution[q] = imapping[s];
	}
	return solution;
}

}

bool use_skew_cyclic(const SparseGraph& g, const KingsGraph& kg){
	assert(kg.width() == kg.height());
	const int n = g.size(), m = kg.height();
	const int vertices_per_path = n / m;
	if(vertices_per_path == 0){ return false; }
	const int length_per_vertex = m / ((n + m - 1) / m);
	if(length_per_vertex == 0){ return false; }
	return true;
}

std::vector<int> skew_cyclic(const SparseGraph& g, const KingsGraph& kg){
	return skew_cyclic_detail::run(g, kg);
}
