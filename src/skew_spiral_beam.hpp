#pragma once
#include <array>
#include <vector>
#include <set>
#include <bitset>
#include <algorithm>
#include <memory>
#include <cstdint>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"

namespace skew_spiral_beam_detail {

static inline bool between(int a, int b, int c){
	return a <= b && b < c;
}

struct Node {
	int score;
	int remaining;
	int vertex;
	int prev;
	std::array<std::shared_ptr<Node>, 13> ancestors;
	std::bitset<512> status;

	Node()
		: score(0)
		, remaining(0)
		, vertex(-1)
		, ancestors()
		, status()
	{ }

	Node(int s, int r, int v, int p, std::shared_ptr<Node> parent)
		: score(s)
		, remaining(r)
		, vertex(v)
		, prev(p)
		, ancestors()
		, status()
	{
		const int n = ancestors.size();
		ancestors[0] = std::move(parent);
		for(int i = 0; ancestors[i] && i + 1 < n; ++i){
			ancestors[i + 1] = ancestors[i]->ancestors[i];
		}
		if(ancestors[0]){ status = ancestors[0]->status; }
		if(v >= 0){ status[v] = true; }
	}

	const Node& kth_ancestor(int k) const {
		const Node *cur = this;
		for(int i = 0; cur && i < ancestors.size(); ++i){
			if(k & (1 << i)){ cur = cur->ancestors[i].get(); }
		}
		return *cur;
	}
};

using NodePtr = std::shared_ptr<Node>;

struct NodeComparator {
	double remaining_weight;
	NodeComparator() : remaining_weight(0.0) { }
	explicit NodeComparator(double rw) : remaining_weight(rw) { }
	bool operator()(const NodePtr& a, const NodePtr& b) const {
		const int sdiff = a->score - b->score;
		const int rdiff = a->remaining - b->remaining;
		return sdiff + rdiff * remaining_weight < 0;
		//if(a->score != b->score){ return a->score < b->score; }
		//return a->remaining < b->remaining;
	}
};

static std::vector<int> make_fill_order(const KingsGraph& kg){
	static const int dy[] = { 0, 1, 0, -1 };
	static const int dx[] = { 1, 0, -1, 0 };
	const int h = kg.height(), w = kg.width();
	std::vector<int> order(h * w, -1);
	for(int i = h * w - 1, y = 0, x = 0, d = 0; i >= 0; --i){
		order[y * w + x] = i;
		const int yy = y + dy[d], xx = x + dx[d];
		if(yy < 0 || h <= yy || xx < 0 || w <= xx || order[yy * w + xx] >= 0){
			d = (d + 1) % 4;
		}
		y += dy[d];
		x += dx[d];
	}
	return order;
}

static std::vector<int> inverse_permutation(const std::vector<int>& p){
	const int n = p.size();
	std::vector<int> q(n);
	for(int i = 0; i < n; ++i){ q[p[i]] = i; }
	return q;
}

static std::vector<int> make_continue_table(const KingsGraph& kg){
	const int h = kg.height(), w = kg.width();
	const auto order = make_fill_order(kg);
	std::vector<int> cont(h * w, -1);
	for(int y = 0; y < h; ++y){
		for(int x = 0; x < w; ++x){
			const int sign = 1 - ((y + x) & 1) * 2;
			const int cur_order = order[y * w + x];
			int prev_order = h * w, prev = -1;
			{	// (+1, +1) or (-1, +1)
				const int yy = y + 1, xx = x + sign;
				if(between(0, yy, h) && between(0, xx, w)){
					const int k = yy * w + xx;
					if(order[k] < cur_order){ prev = k; }
				}
			}
			{	// (-1, -1) or (+1, -1)
				const int yy = y - 1, xx = x - sign;
				if(between(0, yy, h) && between(0, xx, w)){
					const int k = yy * w + xx;
					if(order[k] < cur_order){ prev = k; }
				}
			}
			cont[y * w + x] = prev;
		}
	}
	return cont;
}

static std::vector<uint8_t> make_ignoring_masks(const KingsGraph& kg){
	const int h = kg.height(), w = kg.width();
	const auto cont = make_continue_table(kg);
	std::vector<int> last_table(h * w, -1);
	std::vector<uint8_t> masks(h * w);
	for(int y = 0; y < h; ++y){
		for(int x = 0; x < w; ++x){
			const int q = y * w + x;
			if(cont[q] < 0){ continue; }
			for(const int p : kg[cont[q]]){ last_table[p] = q; }
			const auto neighbors = kg[q];
			for(int i = 0; i < neighbors.size(); ++i){
				if(last_table[neighbors[i]] == q){
					masks[q] |= (1u << i);
				}
			}
		}
	}
	return masks;
}

std::vector<int> run(const SparseGraph& g, const KingsGraph& kg){
	const int beam_width = 100;
	const int n = g.size(), h = kg.height(), w = kg.width();
	const auto inv_order = make_fill_order(kg);
	const auto order = inverse_permutation(inv_order);
	const auto cont = make_continue_table(kg);
	const auto ignore = make_ignoring_masks(kg);

	std::vector<uint8_t> weights(n * n);
	for(int u = 0; u < n; ++u){
		for(const int v : g[u]){ weights[u * n + v] = 1; }
	}

	std::multiset<NodePtr, NodeComparator> states;
	states.insert(std::make_shared<Node>(0, n, -1, -1, nullptr));
	for(int i = 0; i < h * w; ++i){
		decltype(states) next_states(
			NodeComparator(xorshift128() * 2.0 / 0xffffffffu));
		const int q = order[i];
		const auto neighbors = kg[q];
		for(const auto& cptr : states){
			const auto& cur = *cptr;
			const bool force_put = (cur.remaining >= (h * w - i));
			if(!force_put && cont[q] >= 0){
				const int u_step = inv_order[q] - inv_order[cont[q]];
				const int u = cur.kth_ancestor(u_step - 1).vertex;
				if(u >= 0){
					const auto mask = ignore[q];
					int s = cur.score;
					for(size_t j = 0; j < neighbors.size(); ++j){
						if(mask & (1u << j)){ continue; }
						const int p = neighbors[j];
						if(inv_order[p] >= i){ continue; }
						const auto& h = cur.kth_ancestor(i - inv_order[p] - 1);
						s += weights[u * n + h.vertex];
					}
					next_states.insert(
						std::make_shared<Node>(s, cur.remaining, u, 1, cptr));
					if(next_states.size() > beam_width){
						next_states.erase(next_states.begin());
					}
				}
			}
			if(!force_put){
				next_states.insert(std::make_shared<Node>(
					cur.score, cur.remaining, -1, -1, cptr));
				if(next_states.size() > beam_width){
					next_states.erase(next_states.begin());
				}
			}
			for(int u = 0; u < n; ++u){
				if(cur.status[u]){ continue; }
				int s = cur.score;
				for(const int p : kg[q]){
					if(inv_order[p] >= i){ continue; }
					const auto& h = cur.kth_ancestor(i - inv_order[p] - 1);
					const int v = h.vertex;
					if(v >= 0){ s += weights[u * n + v]; }
				}
				next_states.insert(
					std::make_shared<Node>(s, cur.remaining - 1, u, -1, cptr));
				if(next_states.size() > beam_width){
					next_states.erase(next_states.begin());
				}
			}
		}
		states = std::move(next_states);
	}

	std::vector<int> mapping(h * w, -1);
	NodePtr cur = *states.rbegin();
	for(int i = h * w - 1; i >= 0; --i){
		const int v = cur->vertex;
		mapping[order[i]] = v;
		cur = cur->ancestors[0];
	}
	return mapping;
}

}

std::vector<int> skew_spiral_beam(const SparseGraph& g, const KingsGraph& kg){
	return skew_spiral_beam_detail::run(g, kg);
}
