#pragma once

template <typename T>
class PointerRange {
public:
	using value_type = T;
	using const_iterator = const T *;
private:
	const T *m_begin;
	const T *m_end;
public:
	PointerRange(const T *begin, const T *end)
		: m_begin(begin)
		, m_end(end)
	{ }
	bool empty() const { return m_begin == m_end; }
	size_t size() const { return m_end - m_begin; }
	value_type operator[](const size_t i) const { return m_begin[i]; }
	const_iterator cbegin() const { return m_begin; }
	const_iterator cend() const { return m_end; }
	const_iterator begin() const { return m_begin; }
	const_iterator end() const { return m_end; }
};
