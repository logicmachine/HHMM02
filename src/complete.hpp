#pragma once
#include <cassert>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"

bool is_complete_embeddable(const SparseGraph& g, const KingsGraph& kg){
	assert(kg.width() == kg.height());
	return g.size() <= kg.width();
}

std::vector<int> complete_embedding(const SparseGraph& g, const KingsGraph& kg){
	assert(is_complete_embeddable(g, kg));
	const int n = g.size(), w = kg.width(), h = kg.height();
	std::vector<int> mapping(w * h, -1);
	for(int i = 0; i < n; ++i){
		int c = i, d = (i % 2 == 0 ? 1 : -1);
		for(int j = 0; j < n; ++j){
			mapping[j * w + c] = i;
			if(c == 0 && d < 0){
				d = -d;
			}else if(c == n - 1 && d > 0){
				d = -d;
			}else{
				c += d;
			}
		}
	}
	return mapping;
}
