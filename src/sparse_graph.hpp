#pragma once
#include <vector>

class SparseGraph {

private:
	std::vector<std::vector<int>> m_edges;

public:
	SparseGraph()
		: m_edges()
	{ }

	explicit SparseGraph(int n)
		: m_edges(n)
	{ }

	int size() const {
		return m_edges.size();
	}

	void add_edge(int from, int to){
		m_edges[from].push_back(to);
	}

	const std::vector<int>& operator[](int u) const {
		return m_edges[u];
	}

};
