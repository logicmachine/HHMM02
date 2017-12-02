#include <cstdio>
#include <cassert>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <random>
#include <numeric>
#include "../common/rand_gen.h"

void gen_randomgraph(FILE *fp, int V, int E, unsigned long long int seed) {
  Rand rnd_random(seed);
  assert(V - 1 <= E && E <= V * (V-1) / 2);
  fprintf(fp, "%d %d\n", V, E);

  // random device using seed
  std::mt19937_64 engine(seed);
  
  std::vector<int> func(V);
  std::iota(func.begin(), func.end(), 1);
  std::shuffle(func.begin(), func.end(), engine);

  std::set< std::pair<int, int> > edges;
  for(int i=1; i<V; i++) {
    int parent = rnd_random.NextInt(0, i-1);
    int u = func[i], v = func[parent];
    if(u > v) std::swap(u, v);
    fprintf(fp, "%d %d\n", u, v);
    edges.insert( std::make_pair(u, v) );
  }

  for(int i=0; i<E-(V-1); i++) {
    while(1) {
      int u = rnd_random.NextInt(1, V);
      int v = rnd_random.NextInt(1, V);
      if(u == v) continue;
      if(u >  v) std::swap(u, v);

      // do not insert multiple edge
      if(edges.count( std::make_pair(u, v) ) != 0) continue;
      edges.insert( std::make_pair(u, v) );
	
      fprintf(fp, "%d %d\n", u, v);
      break;
    }
  }
}
