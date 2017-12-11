#!/bin/bash
g++ -std=gnu++11 -O2 -o graph_generator graph_generator.cpp
g++ -std=gnu++11 -O2 -o score_evaluator score_evaluator.cpp
for ((i=1000000; i<1000160; i++)); do
	./graph_generator "random${i}.in" "${i}"
done
