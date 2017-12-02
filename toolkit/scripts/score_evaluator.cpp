#include <cstdio>
#include <vector>
#include <queue>
#include <utility>
#include <string>
#include <set>

using Graph  = std::vector< std::vector<int> >;
using Matrix = std::vector< std::vector<int> >;
using Edge   = std::pair<int, int>;

// initial score
const int X = 5000;
// if the edge is preserved, gain A points each
const int A = 100;
// if i-th set size is P, lose B*(P-1) points each
const int B = 1;
// if all edges are preserved, gain C points
const int C = 100000;

void add_edge(Graph &G, int u, int v) {
  G[u].push_back(v);
  G[v].push_back(u);
}

void bfs(const Graph &G, const std::vector<int> &targets, int &cnt_target) {
  int V = (int)G.size();
  std::vector<int> visited(V), is_target(V);
  for(size_t i=0; i<targets.size(); i++) {
    int vertex = targets[i];
    is_target[vertex]++;
  }

  std::queue<int> que;
  int start_point = targets[0];
  que.push(start_point);
  visited[start_point]++;
  cnt_target++;

  while(!que.empty()) {
    int cur_vertex = que.front(); que.pop();
    for(size_t i=0; i<G[cur_vertex].size(); i++) {
      int nxt_vertex = G[cur_vertex][i];
      if(visited[nxt_vertex]) continue;
      if(is_target[nxt_vertex] == 0) continue;

      visited[nxt_vertex]++;
      cnt_target++;
      que.push(nxt_vertex);
    }
  }
}

int main(int argc, char **argv) {
  if(argc != 3) {
    fprintf(stderr, "Usage: %s [input-file] [output-file]\n", argv[0]);
    return 1;
  }

  // read files
  int temp = -1, score = 0;
  bool all_preserved = true;
  bool error_occured = false;
  int add_A = 0, loss_B = 0;
  FILE *in  = fopen(argv[1], "r");
  FILE *out = fopen(argv[2], "r");
  if(in == NULL) {
    fprintf(stderr, "Error (cannot open): %s\n", argv[1]);
    return 1;
  }
  if(out == NULL) {
    fprintf(stderr, "Error (cannot open): %s\n", argv[2]);
    return 1;
  }

  try {
    // read input-file
    // input G
    int V, E;
    fscanf(in, "%d%d", &V, &E);
    Graph G(V);
    std::set<Edge> edges;
    for(int i=0; i<E; i++) {
      int u, v;
      if(fscanf(in, "%d%d", &u, &v) == EOF)   throw 0;
      if(u >= v)                              throw 4;
      u--; v--;
      if(u < 0 || u >= V || v < 0 || v >= V)  throw 2;
      Edge edge = std::make_pair(u, v);
      if(edges.count(edge))                   throw 5;
      edges.insert(edge);
      add_edge(G, u, v);
    }
    
    // input G emb
    int V_emb, E_emb;
    fscanf(in, "%d%d", &V_emb, &E_emb);
    Graph G_emb(V_emb);
    std::set<Edge> edges_emb;
    for(int i=0; i<E_emb; i++) {
      int a, b;
      if(fscanf(in, "%d%d", &a, &b) == EOF)          throw 0;
      if(a >= b)                                     throw 4;
      a--; b--;
      if(a < 0 || a >= V_emb || b < 0 || b >= V_emb) throw 2;
      Edge edge = std::make_pair(a, b);
      if(edges_emb.count(edge))                      throw 5;
      edges_emb.insert(edge);
      add_edge(G_emb, a, b);
    }
    
    // check EOF (input-file)
    if(fscanf(in, "%d", &temp) != EOF)               throw 0;

    // read output-file
    Matrix user_answer(V);
    std::vector<int> used(V_emb, -1);

    for(int i=0; i<V; i++) {
      int num_element;
      if(fscanf(out, "%d", &num_element) == EOF)   throw 1;
      if(num_element < 0 || num_element > V_emb)   throw 3;

      // do not allow empty set
      if(num_element == 0)                         throw 7;

      // the vertex belongs to i-th set
      for(int k=0; k<num_element; k++) {
        int vertex;
        if(fscanf(out, "%d", &vertex) == EOF)      throw 1;
        vertex--;
        if(vertex < 0 || vertex >= V_emb)          throw 3;
        if(used[vertex] >= 0)                      throw 6;
        used[vertex] = i;
        user_answer[i].push_back(vertex);
      }
    }
    
    // check EOF (output-file)
    if(fscanf(out, "%d", &temp) != EOF)            throw 1;
    
    int sum_set_size = 0;
    // connecting check & set-size check
    for(int i=0; i<V; i++) {
      int user_set_size = (int)user_answer[i].size();
      sum_set_size += user_set_size - 1;
      int connect_size = 0;
      bfs(G_emb, user_answer[i], connect_size);
      if(user_set_size != connect_size)            throw 8;
    }

    score = X - B * sum_set_size;
    loss_B = B * sum_set_size;

    // edge preserving check (for all edge \in E)
    for(std::set<Edge>::iterator it=edges.begin(); it!=edges.end(); it++) {
      int a = (*it).first;
      int b = (*it).second;
      
      bool exist_edge = false;
      for(size_t x=0; x<user_answer[a].size(); x++) {
        for(size_t y=0; y<user_answer[b].size(); y++) {
          int u = user_answer[a][x];
          int v = user_answer[b][y];
          if(u >= v) std::swap(u, v);

          Edge cand_edge = std::make_pair(u, v);
          if(edges_emb.count(cand_edge)) {
            exist_edge = true;
            goto judge_edge;
          }
        }
      }
    judge_edge:
      if(exist_edge) {
        score += A;
        add_A += A;
      }
      else {
        all_preserved = false;
      }
    }

    if(all_preserved) score += C;
  }

  // exception
  catch(int mode) {
    error_occured = true;
    // 0, 1
    if(mode < 2) {
      std::string s = (mode ? "output" : "input");
      fprintf(stderr, "Error: the number of lines is invalid (%s)\n", s.c_str());
    }
    // 2, 3
    else if(mode < 4) {
      mode -= 2;
      std::string s = (mode ? "output" : "input");
      fprintf(stderr, "Error: out of range (%s)\n", s.c_str());
    }
    // 4
    else if(mode == 4) {
      fprintf(stderr, "Error: for the edge (u, v), v must be strictly greater than u (input)\n");
    }
    // 5
    else if(mode == 5) {
      fprintf(stderr, "Error: given graph must not have multiple edge (input)\n");
    }
    // 6
    else if(mode == 6) {
      fprintf(stderr, "Error: the number of appearance of vertex x_{i, j} must be at most 1 (output)\n");
    }
    // 7
    else if(mode == 7) {
      fprintf(stderr, "Error: empty set is not allowed (output)\n");
    }
    // 8
    else if(mode == 8) {
      fprintf(stderr, "Error: vertices in set are not connected (output)\n");
    }
  }

  if(error_occured) {
    fprintf(stderr, "There is something wrong with your solution.\n");
  }
  else {
    fprintf(stderr, "Base                 : %7d points\n", X);
    fprintf(stderr, "Add  (preserved edge): %7d points\n", add_A);
    fprintf(stderr, "Loss (size of set)   : %7d points\n", loss_B);
    fprintf(stderr, "Bonus                : %7d points\n\n", (all_preserved ? C : 0));
  }
    fprintf(stderr, "Score                : %7d points\n", score);
  printf("%d\n", score);
  return 0;
}
