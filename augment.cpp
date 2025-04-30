#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/graphsimplify.hpp"
#include "mylib/lwdp.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file]\n", fil);
}

int main(int argc, char **argv){
  if(argc < 4){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  bool reordering = false;
  std::vector<double> pi;
  
  {
    Graph H;
    if(!H.readfromFile(argv[1])){
      fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[1]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    n = H.numV();
    m = H.numE();
    
    std::vector<double> prob(m);
    pi.resize(m);
    {
      FILE *fp;
      if((fp = fopen(argv[2], "r")) == NULL){
        fprintf(stderr, "ERROR: reading probability file %s failed.\n", argv[2]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
      
      for(size_t i=0; i<m; ++i){
        fscanf(fp, "%lf", &prob[i]);
      }
      fclose(fp);
    }
    
    if(argv[3][0] == '!'){
      reordering = true;
      for(size_t i=0; i<m; ++i){
        G.addEdge(H.e[i].first, H.e[i].second);
        pi[i] = prob[i];
      }
    }else{
      if(!G.readfromFile(argv[3])){
        fprintf(stderr, "ERROR: reading order file %s failed.\n", argv[3]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
    
      for(size_t i=0; i<m; ++i){
        pi[i] = prob[H.etovar(G.e[i].first, G.e[i].second)];
      }
    }
  }
  
  
  auto cstart = std::chrono::system_clock::now();
  
  std::vector<AugmentResult> ares;
  LWDP LWSolver;
  LWSolver.solveAugment(G, pi, ares, reordering);
  
  size_t ibest = 0;
  double resbest = 0;
  
  for(size_t i=0; i<m; ++i){
    printf("%3d-%3d : %3d-%3d %.15lf\n", G.e[i].first, G.e[i].second, ares[i].first.first, ares[i].first.second, ares[i].second);
    if(ares[i].second > resbest){
      ibest = i;
      resbest = ares[i].second;
    }
  }
  
  puts("BEST:");
  printf("%3d-%3d : %3d-%3d %.15lf\n", G.e[ibest].first, G.e[ibest].second, ares[ibest].first.first, ares[ibest].first.second, ares[ibest].second);
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}