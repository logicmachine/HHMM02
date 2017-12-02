#pragma once
#include <vector>
#include <cstdint>
#include "pointer_range.hpp"

class KingsGraph {

public:
	using NeighborSet = PointerRange<int16_t>;

private:
	static constexpr int MAX_DEGREE = 8;

	const int m_height;
	const int m_width;

	std::vector<int16_t> m_neighbor_table;

public:
	KingsGraph()
		: m_height(0)
		, m_width(0)
		, m_neighbor_table()
	{ }
	
	KingsGraph(int height, int width)
		: m_height(height)
		, m_width(width)
		, m_neighbor_table(width * height * (MAX_DEGREE + 1))
	{
		static const int dy[] = { -1, -1, -1, 0, 0, 1, 1, 1 };
		static const int dx[] = { -1, 0, 1, -1, 1, -1, 0, 1 };
		for(int uy = 0; uy < height; ++uy){
			for(int ux = 0; ux < width; ++ux){
				const int base = (uy * width + ux) * (MAX_DEGREE + 1);
				for(int i = 0, j = 1; i < 8; ++i){
					const int vy = uy + dy[i], vx = ux + dx[i];
					if(vy < 0 || height <= vy){ continue; }
					if(vx < 0 || width <= vx){ continue; }
					m_neighbor_table[base + j++] = vy * width + vx;
					++m_neighbor_table[base];
				}
			}
		}
	}

	int height() const { return m_height; }
	int width() const { return m_width; }
	int size() const { return m_height * m_width; }

	int row(int v) const { return v / m_width; }
	int col(int v) const { return v % m_width; }

	NeighborSet operator[](int u) const {
		const auto base = m_neighbor_table.data() + u * (MAX_DEGREE + 1);
		return NeighborSet(base + 1, base + 1 + *base);
	}

};
