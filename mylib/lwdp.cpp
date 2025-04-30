#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"
#include "beam_search/beam_search.hpp"
#include "lwdp.hpp"

#include <vector>
#include <utility>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

void LWDP::solve(const Graph& G, const std::vector<double>& pi, std::vector<double>& res, bool reordering){
  allClear();
  
  GraphSimplify GS;
  GS.LWSimplify(G, pi, std::unordered_set<int>());
  
  Graph& Gn = GS.Gnew;
  std::vector<double>& pin = GS.pinew;
  int no = G.numV();
  int mo = G.numE();
  int nn = Gn.numV();
  int mn = Gn.numE();
  std::vector<double> innerres((nn+1)*(nn+1));
  
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
    fprintf(stderr, "Reordering finished.\n");
  }
  
  Gn.buildFrontiers();
  fprintf(stderr, "max_frontier_width = %d\n", Gn.maxFWidth());
  
  innerSolveEmpty(Gn, GS.pinew, innerres);
  
  res.resize((no+1)*(no+1));
  std::fill(res.begin(), res.end(), 0.0);
  res[0] = innerres[0];
  for(int i=1; i<=no; ++i){
    if(GS.oldtonew[i] != 0){
      for(int j=1; j<=no; ++j){
        if(GS.oldtonew[j] != 0) res[i*(no+1)+j] = innerres[GS.oldtonew[i]*(nn+1)+GS.oldtonew[j]];
      }
    }
  }
  
  for(auto itr = GS.hists.rbegin(); itr != GS.hists.rend(); ++itr){
    int x = std::get<0>(*itr);
    int y = std::get<1>(*itr);
    double pri = std::get<2>(*itr);
    for(int i=1; i<=no; ++i){
      res[i*(no+1)+x] = res[x*(no+1)+i] = res[i*(no+1)+y] * pri;
    }
    res[x*(no+1)+x] = 1.0;
  }
}

void LWDP::solveAugment(const Graph& G, const std::vector<double>& pi, std::vector<AugmentResult>& ares, bool reordering){
  std::vector<double> res;
  GraphSimplify GS;
  GS.LWSimplify(G, pi, std::unordered_set<int>());
  
  Graph& Gn = GS.Gnew;
  std::vector<double>& pin = GS.pinew;
  int no = G.numV();
  int mo = G.numE();
  int nn = Gn.numV();
  int mn = Gn.numE();
  std::vector<double> innerres((nn+1)*(nn+1));
  ares.resize(mo);
  
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
    fprintf(stderr, "Reordering finished.\n");
  }
  
  Gn.buildFrontiers();
  fprintf(stderr, "max_frontier_width = %d\n", Gn.maxFWidth());
  
  // Compute the cases when the edge removed by GraphSimlify is augmented
  if(GS.hists.size() > 0){
    allClear();
    innerSolveEmpty(Gn, GS.pinew, innerres);
    
    res.resize((no+1)*(no+1));
    std::fill(res.begin(), res.end(), 0.0);
    res[0] = innerres[0];
    for(int i=1; i<=no; ++i){
      if(GS.oldtonew[i] != 0){
        for(int j=1; j<=no; ++j){
          if(GS.oldtonew[j] != 0) res[i*(no+1)+j] = innerres[GS.oldtonew[i]*(nn+1)+GS.oldtonew[j]];
        }
      }
    }
    
    for(auto itr0 = GS.hists.rbegin(); itr0 != GS.hists.rend(); ++itr0){
      for(auto itr = GS.hists.rbegin(); itr != GS.hists.rend(); ++itr){
        int x = std::get<0>(*itr);
        int y = std::get<1>(*itr);
        double pri = std::get<2>(*itr);
        if(itr0 == itr) pri = 1.0 - (1.0 - pri) * (1.0 - pri);
        for(int i=1; i<=no; ++i){
          res[i*(no+1)+x] = res[x*(no+1)+i] = res[i*(no+1)+y] * pri;
        }
        res[x*(no+1)+x] = 1.0;
      }
      
      int minind = no+3;
      for(int j=1; j<=no; ++j){
        for(int k=j+1; k<=no; ++k){
          if(res[j*(no+1)+k] < res[minind]){
            minind = j*(no+1)+k;
          }
        }
      }
      int evar = G.etovar(std::get<0>(*itr0), std::get<1>(*itr0));
      fprintf(stderr, "<%d>\n", evar);
      ares[evar].first.first = minind / (no+1);
      ares[evar].first.second = minind % (no+1);
      ares[evar].second = res[minind];
    }
  }
  
  // Compute the cases when another edge is augmented
  std::vector<int> newtoold(nn+1);
  for(int i=1; i<=no; ++i){
    if(GS.oldtonew[i] != 0) newtoold[GS.oldtonew[i]] = i;
  }
  
  for(int ii=0; ii<mn; ++ii){
    std::vector<double> piaug(GS.pinew);
    piaug[ii] = 1.0 - (1.0 - piaug[ii]) * (1.0 - piaug[ii]);
    allClear();
    innerSolveEmpty(Gn, piaug, innerres);
    
    res.resize((no+1)*(no+1));
    std::fill(res.begin(), res.end(), 0.0);
    res[0] = innerres[0];
    for(int i=1; i<=no; ++i){
      if(GS.oldtonew[i] != 0){
        for(int j=1; j<=no; ++j){
          if(GS.oldtonew[j] != 0) res[i*(no+1)+j] = innerres[GS.oldtonew[i]*(nn+1)+GS.oldtonew[j]];
        }
      }
    }
    
    for(auto itr = GS.hists.rbegin(); itr != GS.hists.rend(); ++itr){
      int x = std::get<0>(*itr);
      int y = std::get<1>(*itr);
      double pri = std::get<2>(*itr);
      for(int i=1; i<=no; ++i){
        res[i*(no+1)+x] = res[x*(no+1)+i] = res[i*(no+1)+y] * pri;
      }
      res[x*(no+1)+x] = 1.0;
    }
    
    int minind = no+3;
    for(int j=1; j<=no; ++j){
      for(int k=j+1; k<=no; ++k){
        if(res[j*(no+1)+k] < res[minind]){
          minind = j*(no+1)+k;
        }
      }
    }
    int evar = G.etovar(newtoold[Gn.e[ii].first], newtoold[Gn.e[ii].second]);
    fprintf(stderr, "<%d>\n", evar);
    ares[evar].first.first = minind / (no+1);
    ares[evar].first.second = minind % (no+1);
    ares[evar].second = res[minind];
  }
}

void LWDP::innerSolveEmpty(const Graph& G, const std::vector<double>& pi, std::vector<double>& res){
  int n = G.numV();
  int m = G.numE();
  std::vector<int64_t> ssizes(m+1);
  std::vector<size_t> beg(n+1);
  std::vector<size_t> fin(n+1);
  
  maps.resize(m+1);
  dp.resize(m+1);
  res.resize((n+1)*(n+1));
  
  // id=0: root node
  {
    State root;
    root.fill(-1);
    maps[0].emplace(root, 0);
    dp[0].emplace_back(root, 0);
    dp[0][0].p = 1.0;
    ssizes[0] = 1;
  }
  
  for(size_t i=0; i<m; ++i){
    const auto& now_fro = G.fros[i];
    const auto& med_fro = G.mfros[i];
    const auto& next_fro = G.fros[i+1];
    const auto& now_vpos = G.vpos[i];
    const auto& now_ent = G.fent[i];
    const auto& now_lve = G.flve[i];
    size_t kk = now_fro.size();
    size_t tt = med_fro.size();
    size_t ll = next_fro.size();
    dp[i+1].reserve(ssizes[i] * 2);
    maps[i+1].reserve(ssizes[i] * 2);
    ssizes[i+1] = 0;
    
    for(const auto& pos : now_ent) beg[med_fro[pos]] = i;
    for(const auto& pos : now_lve) fin[med_fro[pos]] = i;
    
    for(const auto& ent : maps[i]){
      size_t now_id = ent.second;
      const State& now_state = ent.first;
      
      // generate intermediate state
      State med_state;
      med_state.fill(-1);
      int8_t cc = dp[i][now_id].cnum;
      int8_t cc_old = cc;
      
      memcpy(med_state.data(), now_state.data(), tt);
      for(const auto& pos : now_ent){
        med_state[pos] = cc++;
      }
      
      // lo_state processing
      {
        // generate lo_state
        State lo_state;
        lo_state.fill(-1);
        memcpy(lo_state.data(), med_state.data(), ll);
        for(const auto& pos : now_lve){
          if(pos < ll) lo_state[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(size_t j=0; j<ll; ++j){
          auto& val = lo_state[j];
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
        }
        
        // find or generate id
        size_t lo_id;
        auto it = maps[i+1].find(lo_state);
        if(it != maps[i+1].end()){ // maps[i+1] has already had entry
          lo_id = it->second;
        }else{                     // there is no entry
          maps[i+1].emplace(lo_state, ssizes[i+1]);
          lo_id = ssizes[i+1]++;
          dp[i+1].emplace_back(lo_state, cc_new);
        }
        dp[i][now_id].lo = lo_id;
        
        // vlo equals renum
        std::copy(renum.begin(), renum.begin() + cc_old, dp[i][now_id].vlo.begin());
      }
      
      // hi_state processing
      int8_t cat_to   = med_state[now_vpos.first];
      int8_t cat_from = med_state[now_vpos.second];
      {
        // generate hi_state
        State hi_state;
        hi_state.fill(-1);
        memcpy(hi_state.data(), med_state.data(), ll);
        for(const auto& pos : now_lve){
          if(pos < ll) hi_state[pos] = -1;
        }
        std::vector<int8_t> renum(cc, -1);
        int8_t cc_new = 0;
        for(size_t j=0; j<ll; ++j){
          auto& val = hi_state[j];
          if(val < 0) continue;
          if(renum[val] < 0){
            renum[val] = cc_new++;
            if(val == cat_to)        renum[cat_from] = renum[val];
            else if(val == cat_from) renum[cat_to]   = renum[val];
          }
          val = renum[val];
        }
        
        // find or generate id
        size_t hi_id;
        auto it = maps[i+1].find(hi_state);
        if(it != maps[i+1].end()){ // maps[i+1] has already had entry
          hi_id = it->second;
        }else{                     // there is no entry
          maps[i+1].emplace(hi_state, ssizes[i+1]);
          hi_id = ssizes[i+1]++;
          dp[i+1].emplace_back(hi_state, cc_new);
        }
        dp[i][now_id].hi = hi_id;
        
        // vhi equals renum
        std::copy(renum.begin(), renum.begin() + cc_old, dp[i][now_id].vhi.begin());
      }
    }
    maps[i].clear();
  }
  maps[m].clear();
  
  // dp calculation
  // top-down calc. of p
  for(size_t i=0; i<m; ++i){
    for(const auto& dp_now : dp[i]){
      dp[i+1][dp_now.lo].p += (1.0 - pi[i]) * dp_now.p;
      dp[i+1][dp_now.hi].p += pi[i]         * dp_now.p;
    }
  }
  
  // bottom-up calc. of q
  for(size_t i=m-1; i>=1; --i){
    for(auto&& dp_now : dp[i]){
      for(int8_t a=0; a<dp_now.cnum; ++a){
        for(int8_t b=0; b<a; ++b){
          dp_now.q[a][b] = dp_now.q[b][a];
        }
        dp_now.q[a][a] = 1.0;
        for(int8_t b=a+1; b<dp_now.cnum; ++b){
          if(dp_now.vlo[a] >= 0 && dp_now.vlo[b] >= 0) dp_now.q[a][b]  = (1.0 - pi[i]) * dp[i+1][dp_now.lo].q[dp_now.vlo[a]][dp_now.vlo[b]];
          if(dp_now.vhi[a] >= 0 && dp_now.vhi[b] >= 0) dp_now.q[a][b] += pi[i]         * dp[i+1][dp_now.hi].q[dp_now.vhi[a]][dp_now.vhi[b]];
        }
      }
      if(G.flve[i].size() == 2){
        int8_t cat1 = dp_now.s[G.flve[i][0]];
        int8_t cat2 = dp_now.s[G.flve[i][1]];
        if(cat1 != cat2 && dp_now.vhi[cat1] < 0 && dp_now.vhi[cat2] < 0){
          dp_now.q[cat1][cat2] = dp_now.q[cat2][cat1] = pi[i];
        }
      }      
    }
  }
  
  // sort vertices
  std::vector<int> sortv(n);
  for(int k=0; k<n; ++k) sortv[k] = k+1;
  std::sort(sortv.begin(), sortv.end(), [&fin](const int& a, const int& b){return fin[a] < fin[b];});
  std::vector<size_t> vtopos(n+1);
  for(size_t i=0; i<m; ++i){
    for(auto pos : G.fent[i]){
      vtopos[G.mfros[i][pos]] = pos;
    }
  }
  
  size_t fin0 = fin[sortv[0]];
  size_t pos0 = vtopos[sortv[0]];
  for(size_t k=n-1; k>0; --k){
    int nowv = sortv[k];
    size_t bege = beg[nowv];
    size_t pose = vtopos[nowv];
    
    // bottom-up DP
    if(fin0 <= bege){
      for(auto&& dp_now : dp[bege]){
        for(int8_t a=0; a<dp_now.cnum; ++a){
          dp_now.d[a] = 0.0;
          if(dp_now.vlo[a] >= 0) dp_now.d[a]  = (1.0 - pi[bege]) * dp[bege+1][dp_now.lo].q[dp_now.vlo[a]][dp[bege+1][dp_now.lo].s[pose]];
          if(dp_now.vhi[a] >= 0) dp_now.d[a] += pi[bege]         * dp[bege+1][dp_now.hi].q[dp_now.vhi[a]][dp[bege+1][dp_now.hi].s[pose]];
        }
      }
      for(size_t i=bege-1; i>=fin0; --i){
        for(auto&& dp_now : dp[i]){
          for(int8_t a=0; a<dp_now.cnum; ++a){
            dp_now.d[a] = 0.0;
            if(dp_now.vlo[a] >= 0) dp_now.d[a]  = (1.0 - pi[i]) * dp[i+1][dp_now.lo].d[dp_now.vlo[a]];
            if(dp_now.vhi[a] >= 0) dp_now.d[a] += pi[i]         * dp[i+1][dp_now.hi].d[dp_now.vhi[a]];
          }
        }
      }
    }
    
    // levelwise DP
    res[nowv * (n+1) + nowv] = 1.0;
    for(size_t l=0; l<k; ++l){
      int tmpv = sortv[l];
      int post = vtopos[tmpv];
      size_t shk = std::max(beg[tmpv]+1, fin0);
      for(size_t i=shk+1; i<=fin[tmpv]; ++i){
        if(ssizes[shk] > ssizes[i]) shk = i;
      }
      double tmpres = 0.0;
      if(shk > bege){
        for(const auto& dp_now : dp[shk]){
          tmpres += dp_now.p * dp_now.q[dp_now.s[post]][dp_now.s[pose]];
        }
      }else{
        for(const auto& dp_now : dp[shk]){
          tmpres += dp_now.p * dp_now.d[dp_now.s[post]];
        }
      }
      res[nowv * (n+1) + tmpv] = res[tmpv * (n+1) + nowv] = tmpres;
    }
  }
  res[sortv[0] * (n+1) + sortv[0]] = 1.0;
}