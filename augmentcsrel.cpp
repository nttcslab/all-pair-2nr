#include "mylib/common.hpp"
#include "mylib/graph.hpp"
#include "mylib/graphsimplify.hpp"
#include "mylib/csrel.hpp"
#include "mylib/beam_search/beam_search.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <unordered_set>

using AugmentResult = std::pair<std::pair<int, int>, double>;

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [order_file]\n", fil);
}

int main(int argc, char **argv){
  if(argc < 3){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  bool reordering = false;
  std::vector<double> pi;
  std::unordered_set<int> srcs;
  
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
  
  std::vector<double> res;
  std::vector<AugmentResult> ares(m);
  
  auto cstart = std::chrono::system_clock::now();
  
  GraphSimplify GSLW;
  GSLW.LWSimplify(G, pi, srcs);
  Graph& Gn = GSLW.Gnew;
  std::vector<double>& pin = GSLW.pinew;
  int nn = Gn.numV();
  int mn = Gn.numE();
  
  Graph reoG;
  std::vector<double> reopi;
  
  if(reordering){
    std::vector<Edge> beforder;
    beforder.reserve(mn);
    for(const auto& edg : Gn.e){
      beforder.emplace_back(edg.first-1, edg.second-1);
    }
    std::vector<Edge> aftorder = ordering(nn, beforder);
    for(const auto& edg : aftorder){
      reoG.addEdge(edg.first+1, edg.second+1);
    }
    reopi.resize(mn);
    for(int i=0; i<mn; ++i){
      reopi[i] = pin[Gn.etovar(reoG.e[i].first, reoG.e[i].second)];
    }
    Gn = reoG;
    pin = reopi;
  }
  
  std::vector<double> nrels((nn+1)*(nn+1));
  
  Gn.buildFrontiers();
  fprintf(stderr, "max_frontier_width = %d\n", Gn.maxFWidth());
  //Gn.printFWidth();
  
  CSREL CSSolver;
  
  // Compute the cases when the edge removed by GraphSimlify is augmented
  if(GSLW.hists.size() > 0){
    for(int u=1; u<nn; ++u){
      std::vector<double> csres;
      std::unordered_set<int> tgts = {u};
      
      CSSolver.solve(Gn, pin, tgts, csres);
      
      nrels[u*(nn+1)+u] = 1.0;
      for(int v=u+1; v<=nn; ++v){
        nrels[u*(nn+1)+v] = nrels[v*(nn+1)+u] = csres[v];
      }
    }
    nrels[nn*(nn+1)+nn] = 1.0;
    
    res.resize((n+1)*(n+1));
    std::fill(res.begin(), res.end(), 0.0);
    for(int i=1; i<=n; ++i){
      if(GSLW.oldtonew[i] != 0){
        for(int j=1; j<=n; ++j){
          if(GSLW.oldtonew[j] != 0) res[i*(n+1)+j] = nrels[GSLW.oldtonew[i]*(nn+1)+GSLW.oldtonew[j]];
        }
      }
    }
    
    for(auto itr0 = GSLW.hists.rbegin(); itr0 != GSLW.hists.rend(); ++itr0){
      for(auto itr = GSLW.hists.rbegin(); itr != GSLW.hists.rend(); ++itr){
        int x = std::get<0>(*itr);
        int y = std::get<1>(*itr);
        double pri = std::get<2>(*itr);
        if(itr0 == itr) pri = 1.0 - (1.0 - pri) * (1.0 - pri);
        for(int i=1; i<=n; ++i){
          res[i*(n+1)+x] = res[x*(n+1)+i] = res[i*(n+1)+y] * pri;
        }
        res[x*(n+1)+x] = 1.0;
      }
      
      int minind = n+3;
      for(int j=1; j<=n; ++j){
        for(int k=j+1; k<=n; ++k){
          if(res[j*(n+1)+k] < res[minind]){
            minind = j*(n+1)+k;
          }
        }
      }
      int evar = G.etovar(std::get<0>(*itr0), std::get<1>(*itr0));
      fprintf(stderr, "<%d>\n", evar);
      ares[evar].first.first = minind / (n+1);
      ares[evar].first.second = minind % (n+1);
      ares[evar].second = res[minind];
    }
  }
  
  // Compute the cases when another edge is augmented
  std::vector<int> newtoold(nn+1);
  for(int i=1; i<=n; ++i){
    if(GSLW.oldtonew[i] != 0) newtoold[GSLW.oldtonew[i]] = i;
  }
  
  for(int ii=0; ii<mn; ++ii){
    std::vector<double> piaug(pin);
    piaug[ii] = 1.0 - (1.0 - piaug[ii]) * (1.0 - piaug[ii]);
    
    for(int u=1; u<nn; ++u){
      std::vector<double> csres;
      std::unordered_set<int> tgts = {u};
      
      CSSolver.solve(Gn, piaug, tgts, csres);
      
      nrels[u*(nn+1)+u] = 1.0;
      for(int v=u+1; v<=nn; ++v){
        nrels[u*(nn+1)+v] = nrels[v*(nn+1)+u] = csres[v];
      }
    }
    nrels[nn*(nn+1)+nn] = 1.0;
    
    res.resize((n+1)*(n+1));
    std::fill(res.begin(), res.end(), 0.0);
    for(int i=1; i<=n; ++i){
      if(GSLW.oldtonew[i] != 0){
        for(int j=1; j<=n; ++j){
          if(GSLW.oldtonew[j] != 0) res[i*(n+1)+j] = nrels[GSLW.oldtonew[i]*(nn+1)+GSLW.oldtonew[j]];
        }
      }
    }
    
    for(auto itr = GSLW.hists.rbegin(); itr != GSLW.hists.rend(); ++itr){
      int x = std::get<0>(*itr);
      int y = std::get<1>(*itr);
      double pri = std::get<2>(*itr);
      for(int i=1; i<=n; ++i){
        res[i*(n+1)+x] = res[x*(n+1)+i] = res[i*(n+1)+y] * pri;
      }
      res[x*(n+1)+x] = 1.0;
    }
    
    int minind = n+3;
    for(int j=1; j<=n; ++j){
      for(int k=j+1; k<=n; ++k){
        if(res[j*(n+1)+k] < res[minind]){
          minind = j*(n+1)+k;
        }
      }
    }
    int evar = G.etovar(newtoold[Gn.e[ii].first], newtoold[Gn.e[ii].second]);
    fprintf(stderr, "<%d>\n", evar);
    ares[evar].first.first = minind / (n+1);
    ares[evar].first.second = minind % (n+1);
    ares[evar].second = res[minind];
  }
  
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