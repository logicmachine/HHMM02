#pragma once
#include <vector>
#include <algorithm>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"
#include "xorshift128.hpp"

namespace skew_annealing_detail {

class State {

private:
	const SparseGraph& m_g;
	const KingsGraph& m_kg;

	std::vector<int> m_mapping;
	std::vector<int> m_imapping;

public:

};

}

std::vector<int> skew_annealing(const SparseGraph& g, const KingsGraph& kg){
	return mapping;
}
