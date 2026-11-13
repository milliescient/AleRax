#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <vector>

/**
 *  Ascertainment correction for the species coverage filter of the undated
 *  DL model: the likelihood of a family is conditioned on the family being
 *  observed in at least k species (k = --min-covered-species), the filter
 *  that trimFamilies applies, instead of on survival alone.
 *
 *  Let E^C_e be the probability that a gene on species branch e leaves
 *  observed copies only in species of the leaf set C. It is the smaller root
 *  of E = pL + pS x + pD E^2, where x = E^C_f E^C_g on an internal branch, and
 *  x = 1 (s in C) or fm_s (s not in C) on a terminal branch s. With
 *  P(C) = sum_e OP_e E^C_e, inclusion-exclusion gives
 *    P(#species < k) = sum_{|C|<k} P(C) sum_{j=0}^{k-1-|C|} (-1)^j C(n-|C|,j),
 *  and the factor is 1 - P(#species < k). C = {} is the usual extinction.
 */
namespace CoverageConditioning {

/**
 *  Minimum number of covered species to condition on; 0 means off
 *  (condition on survival, as upstream)
 */
inline unsigned int &minSpecies() {
  static unsigned int k = 0;
  return k;
}

inline bool enabled() { return minSpecies() > 0; }

/**
 *  Smaller root of E = a + pD E^2 with a = pL + pS x, in a form that is
 *  stable as pD goes to 0
 */
inline double extinctionRoot(double pD, double a) {
  double disc = 1.0 - 4.0 * pD * a;
  return 2.0 * a / (1.0 + std::sqrt(disc > 0.0 ? disc : 0.0));
}

/**
 *  The species tree in postorder (children before parents): left and right
 *  child positions in that order (-1 for a leaf), the event probabilities,
 *  origination probabilities and fractions missing at each position.
 */
struct Tree {
  std::vector<int> left, right;
  std::vector<double> pD, pL, pS, op, fm;
};

/**
 *  Extinction with the leaves flagged in `inC` unconstrained
 */
inline void maskedExtinction(const Tree &t, const std::vector<char> &inC,
                             std::vector<double> &E) {
  auto n = t.left.size();
  E.resize(n);
  for (size_t i = 0; i < n; ++i) {
    double x;
    if (t.left[i] < 0) {
      x = inC[i] ? 1.0 : t.fm[i];
    } else {
      x = E[t.left[i]] * E[t.right[i]];
    }
    E[i] = extinctionRoot(t.pD[i], t.pL[i] + t.pS[i] * x);
  }
}

inline double binom(unsigned int n, unsigned int j) {
  if (j > n) {
    return 0.0;
  }
  double r = 1.0;
  for (unsigned int i = 1; i <= j; ++i) {
    r = r * static_cast<double>(n - j + i) / static_cast<double>(i);
  }
  return r;
}

/**
 *  Probability that a family is observed in at least k species
 */
inline double coverageProbability(const Tree &t, unsigned int k) {
  std::vector<size_t> leaves;
  for (size_t i = 0; i < t.left.size(); ++i) {
    if (t.left[i] < 0) {
      leaves.push_back(i);
    }
  }
  auto n = static_cast<unsigned int>(leaves.size());
  if (k > n) {
    return 0.0;
  }
  std::vector<double> coef(k, 0.0);
  for (unsigned int c = 0; c < k; ++c) {
    for (unsigned int j = 0; j + c < k; ++j) {
      coef[c] += ((j % 2) ? -1.0 : 1.0) * binom(n - c, j);
    }
  }
  std::vector<char> inC(t.left.size(), 0);
  std::vector<double> E;
  double below = 0.0; // P(#species < k)
  // every leaf subset of size < k, each visited once
  std::function<void(size_t, unsigned int)> visit = [&](size_t start,
                                                        unsigned int size) {
    maskedExtinction(t, inC, E);
    double p = 0.0;
    for (size_t i = 0; i < E.size(); ++i) {
      p += t.op[i] * E[i];
    }
    below += coef[size] * p;
    if (size + 1 >= k) {
      return;
    }
    for (size_t i = start; i < n; ++i) {
      inC[leaves[i]] = 1;
      visit(i + 1, size + 1);
      inC[leaves[i]] = 0;
    }
  };
  visit(0, 0);
  double pass = 1.0 - below;
  return pass > 0.0 ? pass : 0.0;
}

/**
 *  coverageProbability with a cache shared by all families of this rank:
 *  under global or per-species rates every family has the same factor
 */
inline double cachedCoverageProbability(const Tree &t, unsigned int k) {
  static std::map<std::vector<double>, double> cache;
  std::vector<double> key;
  key.reserve(7 * t.left.size() + 1);
  key.push_back(static_cast<double>(k));
  for (size_t i = 0; i < t.left.size(); ++i) {
    key.push_back(t.left[i]);
    key.push_back(t.right[i]);
    key.push_back(t.pD[i]);
    key.push_back(t.pL[i]);
    key.push_back(t.pS[i]);
    key.push_back(t.op[i]);
    key.push_back(t.fm[i]);
  }
  auto it = cache.find(key);
  if (it != cache.end()) {
    return it->second;
  }
  if (cache.size() > 100000) {
    cache.clear();
  }
  double p = coverageProbability(t, k);
  cache.emplace(std::move(key), p);
  return p;
}

} // namespace CoverageConditioning
