#pragma once

template <typename T>
class Matrix {

public:
	using value_type = T;

private:
	size_t m_rows;
	size_t m_columns;
	std::vector<value_type> m_data;

public:
	Matrix()
		: m_rows(0)
		, m_columns(0)
		, m_data()
	{ }

	Matrix(size_t r, size_t c, const value_type& x = value_type())
		: m_rows(r)
		, m_columns(c)
		, m_data(r * c, x)
	{ }


	size_t rows() const {
		return m_rows;
	}

	size_t columns() const {
		return m_columns;
	}


	const value_type& operator()(size_t i, size_t j) const {
		return m_data[i * m_columns + j];
	}

	value_type& operator()(size_t i, size_t j){
		return m_data[i * m_columns + j];
	}

};

