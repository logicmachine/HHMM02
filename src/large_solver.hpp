#pragma once
#include <vector>
#include <algorithm>
#include <memory>
#include <bitset>
#include <chrono>
#include <cstdint>
#include <cassert>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"
#include "matrix.hpp"
#include "pointer_range.hpp"
#include "radix_sort.hpp"
#include "xorshift128.hpp"

#define DEBUG_DUMP

namespace large_solver_detail {

//---------------------------------------------------------------------------
// Utilities
//---------------------------------------------------------------------------
template <typename T>
static std::vector<T> inverse_permutation(const std::vector<T>& p){
	const size_t n = p.size();
	std::vector<T> q(n);
	for(size_t i = 0; i < n; ++i){ q[p[i]] = i; }
	return q;
}

template <typename T>
static Matrix<T> make_dense_graph(const SparseGraph& g){
	const int n = g.size();
	Matrix<T> mat(n, n);
	for(int u = 0; u < n; ++u){
		for(const int v : g[u]){ mat(u, v) = 1; }
	}
	return mat;
}

class EmbedGraph {

private:
	int m_size;
	std::vector<int> m_rowptr;
	std::vector<int> m_colind;
	std::vector<int> m_mapping;

	std::vector<int> make_next_table(const KingsGraph& kg){
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

	std::vector<int> make_segment_mapping(int n, const KingsGraph& kg){
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
				if(n % m == 0 || i < n % m || j != (n_ranges - 1) / 2){
					++v;
				}
			}
		}
		return mapping;
	}

	SparseGraph make_segment_adjacency_list(
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

	std::vector<int> compute_placement_order(
		const SparseGraph& g_emb,
		const std::vector<int>& mapping)
	{
		const int n = g_emb.size();
		const int m = static_cast<int>(sqrt(mapping.size()));
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
			const int q =
				max_element(degrees.begin(), degrees.end()) - degrees.begin();
			order.push_back(q);
			degrees[q] = -1;
			for(const int p : g_emb[q]){
				if(degrees[p] >= 0){ ++degrees[p]; }
			}
		}
		return order;
	}

public:
	EmbedGraph()
		: m_size(0)
		, m_rowptr(1)
		, m_colind()
		, m_mapping()
	{ }

	EmbedGraph(int n, const KingsGraph& kg)
		: m_size(n)
		, m_rowptr(n + 1)
		, m_colind()
		, m_mapping(kg.size())
	{
		assert(kg.height() == kg.width());
		const auto mapping = make_segment_mapping(n, kg);
		const auto adjacency = make_segment_adjacency_list(mapping, kg);
		const auto order = compute_placement_order(adjacency, mapping);
		const auto inv_order = inverse_permutation(order);
		for(int i = 0; i < n; ++i){
			const int u = order[i];
			for(const auto& v : adjacency[u]){
				m_colind.push_back(inv_order[v]);
			}
			m_rowptr[i + 1] = static_cast<int>(m_colind.size());
		}
		for(int i = 0; i < kg.size(); ++i){
			m_mapping[i] = inv_order[mapping[i]];
		}
	}

	int size() const {
		return m_size;
	}

	PointerRange<int> operator[](int i) const {
		return PointerRange<int>(
			m_colind.data() + m_rowptr[i],
			m_colind.data() + m_rowptr[i + 1]);
	}

	int mapping(int v) const {
		return m_mapping[v];
	}

};


//---------------------------------------------------------------------------
// Beam search
//---------------------------------------------------------------------------
template <size_t MAX_N>
class BeamSearch {

private:
	struct Node {
		int score;
		std::bitset<MAX_N> status;
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

public:
	std::vector<int> operator()(
		const SparseGraph& g, const EmbedGraph& eg, size_t beam_width)
	{
		assert(g.size() == eg.size());
		const int n = g.size();
		const auto g_mat = make_dense_graph<uint8_t>(g);

		RadixSortEngine<int, std::pair<int, int>> radix_sort;
		std::vector<Node> states(1);
		for(int q = 0; q < n; ++q){
			std::vector<int> neighbors;
			for(const auto p : eg[q]){
				if(p < q){ neighbors.push_back(p); }
			}

			std::vector<int> next_scores;
			std::vector<std::pair<int, int>> next_moves;
			next_scores.reserve(states.size() * (n - q));
			next_moves.reserve(states.size() * (n - q));

			for(int i = 0; i < states.size(); ++i){
				const auto& cur = states[i];
				for(int u = 0; u < n; ++u){
					if(cur.status[u]){ continue; }
					int s = cur.score;
					for(const int p : neighbors){
						s += g_mat(u, cur.history[p]);
					}
					next_scores.push_back(s);
					next_moves.emplace_back(i, u);
				}
			}

			const size_t n_next = next_scores.size();
			radix_sort(n_next, next_scores.data(), next_moves.data());
			const size_t first = n_next - std::min(n_next, beam_width);
			std::vector<Node> next_states;
			next_states.reserve(n_next - first);
			for(size_t i = first; i < n_next; ++i){
				const int score = next_scores[i];
				const int j = next_moves[i].first;
				const int u = next_moves[i].second;
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

		const auto& history = states.back().history;
		std::vector<int> mapping(n);
		for(int i = 0; i < n; ++i){
			mapping[history[i]] = i;
		}
#ifdef DEBUG_DUMP
		std::cerr << "BS: " << states.back().score << std::endl;
#endif
		return mapping;
	}

};


//---------------------------------------------------------------------------
// Simulated Annealing
//---------------------------------------------------------------------------
template <size_t MAX_N>
class SimulatedAnnealing {

public:
	std::vector<int> operator()(
		const SparseGraph& g,
		const EmbedGraph& eg,
		const std::vector<int>& init_mapping,
		const std::chrono::duration<double> time_limit)
	{
		const int n = g.size();
		const auto g_mat = make_dense_graph<uint8_t>(g);

		int score = 0;
		auto mapping = init_mapping;
		auto imapping = inverse_permutation(mapping);

		int best_score = 0;
		auto best_mapping = mapping;

		const auto start_time = std::chrono::steady_clock::now();
		auto last_time = start_time;
		auto max_duration = last_time - start_time;

		const double max_temp = 0.3;
		const double min_temp = 0.1;
		do {
			const double progress = (last_time - start_time) / time_limit;
			const double temperature =
				max_temp - (max_temp - min_temp) * progress;

			const double e = exp(-1.0 / temperature);
			std::array<uint32_t, 256> threshold_table;
			double exp_acc = 1.0;
			for(int i = 0; i < threshold_table.size(); ++i){
				threshold_table[i] =
					static_cast<uint32_t>(exp_acc * 0xffffffffu);
				exp_acc *= e;
			}

			for(int iter = 0; iter < 1000; ++iter){
				const int w = modulus_random(n);
				for(const int v : g[w]){
					const auto neighbors = eg[mapping[w]];
					const int q =
						neighbors[modulus_random(neighbors.size())];
					const int p = mapping[v];
					if(p == q){ continue; }
					const int u = imapping[q];
					int diff = 0;
					for(const int r : eg[q]){
						if(r != p){ diff -= g_mat(u, imapping[r]); }
					}
					for(const int r : eg[p]){
						diff += g_mat(u, imapping[r]);
					}
					for(const int r : eg[p]){
						if(r != q){ diff -= g_mat(v, imapping[r]); }
					}
					for(const int r : eg[q]){
						diff += g_mat(v, imapping[r]);
					}
					const bool update_flag =
						threshold_table[std::max(0, -diff)] >= xorshift128();
					if(update_flag){
						mapping[u] = p;
						mapping[v] = q;
						imapping[q] = v;
						imapping[p] = u;
						score += diff;
						if(score > best_score){
							// std::cerr << "SA: +" << (score - best_score) << " (acc: +" << score << ")" << std::endl;
							best_score = score;
							best_mapping = mapping;
						}
					}
				}
			}

			const auto now = std::chrono::steady_clock::now();
			max_duration = std::max(max_duration, now - last_time);
			last_time = now;
		} while((last_time - start_time) + max_duration * 1.2 < time_limit);

#ifdef DEBUG_DUMP
		std::cerr << "SA: +" << best_score << std::endl;
#endif
		return best_mapping;
	}

};


//---------------------------------------------------------------------------
// Hill climbing
//---------------------------------------------------------------------------
class HillClimbing {

public:
	std::vector<int> operator()(
		const SparseGraph& g,
		const EmbedGraph& eg,
		const std::vector<int>& init_mapping)
	{
		const int n = g.size();
		const auto g_mat = make_dense_graph<uint8_t>(g);
		auto imapping = inverse_permutation(init_mapping);
#ifdef DEBUG_DUMP
		int improve = 0;
#endif
		while(true){
			int best_q = -1, best_p = -1, best_score = 0;
			for(int q = 0; q < n; ++q){
				const int u = imapping[q];
				for(int p = q + 1; p < n; ++p){
					const int v = imapping[p];
					int score = 0;
					for(const int r : eg[q]){
						if(r != p){ score -= g_mat(u, imapping[r]); }
					}
					for(const int r : eg[p]){
						score += g_mat(u, imapping[r]);
					}
					for(const int r : eg[p]){
						if(r != q){ score -= g_mat(v, imapping[r]); }
					}
					for(const int r : eg[q]){
						score += g_mat(v, imapping[r]);
					}
					if(score > best_score){
						best_score = score;
						best_q = q;
						best_p = p;
					}
				}
			}
			if(best_score <= 0){ break; }
			std::swap(imapping[best_q], imapping[best_p]);
#ifdef DEBUG_DUMP
			improve += best_score;
#endif
		}

#ifdef DEBUG_DUMP
		std::cerr << "HC: +" << improve << std::endl;
#endif
		return inverse_permutation(imapping);
	}

};


//---------------------------------------------------------------------------
// Solver
//---------------------------------------------------------------------------
template <size_t MAX_N>
std::vector<int> run(const SparseGraph& g, const KingsGraph& kg){
	assert(kg.width() == kg.height());
	EmbedGraph eg(g.size(), kg);

	BeamSearch<MAX_N> beam_search;
	std::vector<int> mapping = beam_search(g, eg, 4000);

	SimulatedAnnealing<MAX_N> simulated_annealing;
	mapping = simulated_annealing(g, eg, mapping, std::chrono::milliseconds(10000));

	HillClimbing hill_climbing;
	mapping = hill_climbing(g, eg, mapping);

	const auto imapping = inverse_permutation(mapping);
	std::vector<int> solution(kg.size(), -1);
	for(int q = 0; q < kg.size(); ++q){
		if(eg.mapping(q) < 0){ continue; }
		solution[q] = imapping[eg.mapping(q)];
	}
	return solution;
}

}

std::vector<int> solve_large(const SparseGraph& g, const KingsGraph& kg){
	const int n = g.size();
	if(n <= 64 * 1){ return large_solver_detail::run<64 * 1>(g, kg); }
	if(n <= 64 * 2){ return large_solver_detail::run<64 * 2>(g, kg); }
	if(n <= 64 * 3){ return large_solver_detail::run<64 * 3>(g, kg); }
	if(n <= 64 * 4){ return large_solver_detail::run<64 * 4>(g, kg); }
	if(n <= 64 * 5){ return large_solver_detail::run<64 * 5>(g, kg); }
	if(n <= 64 * 6){ return large_solver_detail::run<64 * 6>(g, kg); }
	if(n <= 64 * 7){ return large_solver_detail::run<64 * 7>(g, kg); }
	if(n <= 64 * 8){ return large_solver_detail::run<64 * 8>(g, kg); }
	return std::vector<int>();
}
