#pragma once
#include <array>
#include <vector>
#include <utility>

template <typename Key, typename Value>
class RadixSortEngine {

public:
	using key_type = Key;
	using value_type = Value;

private:
	static const size_t step_size = 256u;
	static const size_t bucket_size = 16u;
	using key_bucket_type = std::array<key_type, bucket_size>;
	using value_bucket_type = std::array<value_type, bucket_size>;

	std::array<size_t, step_size> m_histogram;
	std::array<size_t, step_size> m_counters;
	std::vector<std::pair<key_bucket_type, value_bucket_type>> m_buckets;
	std::vector<key_type> m_key_buffer;
	std::vector<value_type> m_value_buffer;

public:
	RadixSortEngine()
		: m_buckets(step_size)
	{ }

	void operator()(size_t n, key_type *key, value_type *values){
		if(n > m_key_buffer.size()){ m_key_buffer.resize(n); }
		if(n > m_value_buffer.size()){ m_value_buffer.resize(n); }
		key_type *keys_src = key, *keys_dest = m_key_buffer.data();
		value_type *values_src = values, *values_dest = m_value_buffer.data();
		for(size_t d = 0; d < sizeof(key_type); ++d){
			const size_t shift = d * 8u;
			for(size_t i = 0; i < step_size; ++i){
				m_histogram[i] = 0u;
				m_counters[i] = 0u;
			}
			for(size_t i = 0; i < n; ++i){
				const size_t k = (keys_src[i] >> shift) & (step_size - 1);
				++m_histogram[k];
			}
			for(size_t i = 0, s = 0; i < step_size; ++i){
				const size_t t = m_histogram[i];
				m_histogram[i] = s;
				s += t;
			}
			for(size_t i = 0; i < n; ++i){
				const size_t k = (keys_src[i] >> shift) & (step_size - 1);
				const size_t p = m_counters[k]++;
				m_buckets[k].first[p] = std::move(keys_src[i]);
				m_buckets[k].second[p] = std::move(values_src[i]);
				if(m_counters[k] == bucket_size){
					const size_t base = m_histogram[k];
					for(size_t j = 0; j < m_counters[k]; ++j){
						keys_dest[base + j] = std::move(m_buckets[k].first[j]);
						values_dest[base + j] = std::move(m_buckets[k].second[j]);
					}
					m_histogram[k] += m_counters[k];
					m_counters[k] = 0;
				}
			}
			for(size_t k = 0; k < step_size; ++k){
				const size_t base = m_histogram[k];
				for(size_t j = 0; j < m_counters[k]; ++j){
					keys_dest[base + j] = std::move(m_buckets[k].first[j]);
					values_dest[base + j] = std::move(m_buckets[k].second[j]);
				}
			}
			std::swap(keys_src, keys_dest);
			std::swap(values_src, values_dest);
		}
		if(sizeof(key_type) % 2 == 1){
			for(size_t i = 0; i < n; ++i){
				keys_dest[i] = std::move(keys_src[i]);
				values_dest[i] = std::move(values_src[i]);
			}
		}
	}

};
