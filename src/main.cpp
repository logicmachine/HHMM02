#include <iostream>
#include <cmath>
#include "sparse_graph.hpp"
#include "kings_graph.hpp"
#include "complete.hpp"
#include "spiral_greedy.hpp"

int main(){
	std::ios_base::sync_with_stdio(false);
	int n, m;
	std::cin >> n >> m;
	SparseGraph g(n);
	for(int i = 0; i < m; ++i){
		int u, v;
		std::cin >> u >> v;
		--u; --v;
		g.add_edge(u, v);
		g.add_edge(v, u);
	}
	int n_emb, m_emb;
	std::cin >> n_emb >> m_emb;
	for(int i = 0; i < m_emb; ++i){
		int u, v;
		std::cin >> u >> v;
	}
	const int wh_emb = static_cast<int>(sqrt(n_emb));
	KingsGraph kg(wh_emb, wh_emb);

	std::vector<int> mapping(n_emb, -1);
	if(is_complete_embeddable(g, kg)){
		mapping = complete_embedding(g, kg);
	}else{
		mapping = spiral_greedy(g, kg);
	}

	std::vector<std::vector<int>> imapping(n);
	for(int i = 0; i < n_emb; ++i){
		if(mapping[i] >= 0){ imapping[mapping[i]].push_back(i); }
	}
	for(int i = 0; i < n; ++i){
		std::cout << (i + 1);
		for(const int x : imapping[i]){
			std::cout << " " << (x + 1);
		}
		std::cout << std::endl;
	}
	return 0;
}
