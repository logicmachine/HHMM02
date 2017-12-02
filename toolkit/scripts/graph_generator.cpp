#include <cstdio>
#include <cstring>
#include <algorithm>
#include <random>
#include "../common/rand_gen.h"
#include "../common/constraints.h"
#include "../graph/kingsgraph_square.h"
#include "../graph/randomgraph.h"

unsigned long long int str2ull(char *str) {
  char* end_point;
  unsigned long long int ret_number = strtoull(str, &end_point, 10);
  if(*end_point != '\0') {
    fprintf(stderr, "Error: invalid input: \"%s\"\n", str);
    exit(1);
  }
  return ret_number;
}

int main(int argc, char** argv) {
  if(argc != 3) {
    fprintf(stderr, "Usage: %s [filename] [seed]\n", argv[0]);
    return 1;
  }

  FILE *fp = fopen(argv[1], "w");
  if(fp == NULL) {
    fprintf(stderr, "Error: cannot open the file\n");
    return 1;
  }

  unsigned long long int seed = str2ull(argv[2]);
  Rand rnd(seed);
  int V = rnd.NextInt(MIN_V, MAX_V);
  int E = rnd.NextInt(V-1, std::min(V*(V-1)/2, MAX_E));
  gen_randomgraph(fp, V, E, seed);
      
  int N;
  while(1) {
    N = rnd.NextInt(MIN_N, MAX_N);
    if(N*N >= V) break;
  }
  gen_kingsgraph(fp, N);
  fclose(fp);
  return 0;
}
