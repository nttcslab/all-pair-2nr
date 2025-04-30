#ifndef REL_LWDP_HPP
#define REL_LWDP_HPP

#include "common.hpp"
#include "graph.hpp"
#include "graphsimplify.hpp"

#include <array>
#include <vector>
#include <utility>
#include <cassert>
#include <tuple>
#include <unordered_set>
#include <unordered_map>

using State = std::array<int8_t, 16>;

namespace std{
  template <>
  struct hash<State>{
    public:
    uint64_t operator()(const State& s) const{
      const uint64_t *p = reinterpret_cast<const uint64_t*>(s.data());
      return *p * 314159257ULL + *(p+1);
    }
  };
}

class DPBlock{
public:
  State s;
  double p = 0.0;
  std::vector<std::vector<double>> q;
  std::vector<double> d;
  int8_t cnum;
  size_t lo;
  size_t hi;
  std::vector<int8_t> vlo;
  std::vector<int8_t> vhi;
  
  DPBlock() {};
  DPBlock(const State& _s, int8_t _siz): s(_s), cnum(_siz), q(_siz), d(_siz), vlo(_siz), vhi(_siz) {
    for(auto&& v : q) v.resize(_siz);
  };
  DPBlock(State&& _s, int8_t _siz): s(_s), cnum(_siz), q(_siz), d(_siz), vlo(_siz), vhi(_siz) {
    for(auto&& v : q) v.resize(_siz);
  };
};

using AugmentResult = std::pair<std::pair<int, int>, double>;

class LWDP{
public:
  std::vector<std::vector<DPBlock>> dp;
  std::vector<std::unordered_map<State, size_t>> maps;
  std::vector<std::pair<int, int8_t>> vtoipos;
  
  LWDP(){};
  
  void solve(const Graph& G, const std::vector<double>& pi, std::vector<double>& res, bool reordering);
  void solveAugment(const Graph& G, const std::vector<double>& pi, std::vector<AugmentResult>& res, bool reordering);
  void innerSolveEmpty(const Graph& G, const std::vector<double>& pi, std::vector<double>& res);
  void allClear(){
    for(auto&& ret : dp) ret.clear();
    dp.clear();
  }
};

#endif // REL_LWDP_HPP
