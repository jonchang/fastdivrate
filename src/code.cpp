#include <Rcpp.h>
using namespace Rcpp;

std::size_t c_which(NumericVector haystack, std::size_t needle){
  // assumes that the needle is in fact in the haystack
  int nx = haystack.size();
  for(int i = 0; i < nx; i++) {
    if (haystack[i] == needle) {
      return i;
    }
  }
  return 0;
}

//' DR statistic (Rcpp version)
//'
//' Computes the Jetz DR rate
//'
//' @param x an ape::phylo object
//'
//' @return vector of rates for each tip
//' @export
//' @examples
//' tree <- ape::rcoal(1000)
//' DR_statistic_C(tree)
// [[Rcpp::export]]
NumericVector DR_statistic_C(List x) {
  if (!x.inherits("phylo")) stop("Input must be a phylogenetic tree");

  CharacterVector tips = as<CharacterVector>(x["tip.label"]);
  std::size_t ntips = tips.size();
  std::size_t rootnode = ntips + 1;
  NumericVector sprates(ntips, 0.0);
  sprates.names() = tips;
  // save a copy of these for constant-time access instead of
  // multiple pointer derefs in a hot loop
  NumericVector edge_length = as<NumericVector>(x["edge.length"]);
  NumericMatrix edge = as<NumericMatrix>(x["edge"]);
  NumericVector edge1 = edge(_, 0);
  NumericVector edge2 = edge(_, 1);
  for (std::size_t ii = 0; ii != ntips; ++ii) {
    // Although Rcpp is 0 offset, edges in the phylo edge matrix use R numbering
    std::size_t node = ii + 1;
    std::size_t index = 1;
    double qx = 0;
    while (node != rootnode) {
      // avoid horrible incompatible template errors by extracting a scalar value
      // instead of relying on magical Rcpp sugared subsets
      std::size_t parent = c_which(edge2, node);
      double el = edge_length[parent];
      node = edge1[parent];
      qx = qx + el * (1 / pow(2, index - 1));
      ++index;
    }
    sprates[ii] = 1 / qx;
  }
  return sprates;
}

NumericVector stl_sort(NumericVector x) {
    NumericVector y = clone(x);
    std::sort(y.begin(), y.end());
    return y;
}


// TODO: make the R function do fancy dispatching based on the arguments passed
//' DR statistic with sampling times
//'
//' Computes the Jetz DR rate, but sample at certain times in the psat
//'
//' @param x an ape::phylo object
//' @param sample_times a vector of ages
//'
//' @return matrix of rates for each tip, by sample time
//' @export
// [[Rcpp::export]]
NumericMatrix DR_statistic_time(List x, NumericVector sample_times) {
    if (!x.inherits("phylo")) stop("Input must be a phylogenetic tree");
    sample_times = stl_sort(sample_times);

    CharacterVector tips = as<CharacterVector>(x["tip.label"]);
    std::size_t ntips = tips.size();
    std::size_t rootnode = ntips + 1;

    NumericMatrix allrates(ntips, sample_times.size());
    rownames(allrates) = tips;
    colnames(allrates) = as<CharacterVector>(sample_times);

    // save a copy of these for constant-time access instead of
    // multiple pointer derefs in a hot loop
    NumericVector edge_length = as<NumericVector>(x["edge.length"]);
    NumericMatrix edge = as<NumericMatrix>(x["edge"]);
    NumericVector edge1 = edge(_, 0);
    NumericVector edge2 = edge(_, 1);

    for (std::size_t ii = 0; ii != ntips; ++ii) {
        // Although Rcpp is 0 offset, edges in the phylo edge matrix use R numbering
        std::size_t node = ii + 1;
        std::size_t index = 1;
        double cur_time = 0;
        std::size_t sample_idx = 0;
        double qx = 0;
        while (node != rootnode) {
            // avoid horrible incompatible template errors by extracting a scalar value
            // instead of relying on magical Rcpp sugared subsets
            std::size_t parent = c_which(edge2, node);
            double el = edge_length[parent];
            cur_time += el;
            if (sample_times[sample_idx] < cur_time) {
                allrates(ii, sample_idx) = 1 / qx;
                sample_idx++;
            }
            qx = qx + el * (1 / pow(2, index - 1));
            ++index;
            node = edge1[parent];
            if (node == rootnode) {
                // add in the last count if we've reached the end of the tree but still have
                // a sampled time that hasn't been included
                while (sample_idx < sample_times.size()) {
                    allrates(ii, sample_idx) = 1 / qx;
                    sample_idx++;
                }
            }
        }
    }
    return allrates;
}

//' Interval node statistic, with sampling times
//'
//' Computes the interval node statistic, sampling at certain times in the psat
//'
//' @param x an ape::phylo object
//' @param sample_times a vector of ages
//'
//' @return matrix of rates for each tip, by sample time
//' @export
// [[Rcpp::export]]
NumericMatrix interval_node_statistic(List x, NumericVector sample_times) {
    if (!x.inherits("phylo")) stop("Input must be a phylogenetic tree");
    sample_times = stl_sort(sample_times);

    CharacterVector tips = as<CharacterVector>(x["tip.label"]);
    std::size_t ntips = tips.size();
    std::size_t rootnode = ntips + 1;

    NumericMatrix allrates(ntips, sample_times.size());
    rownames(allrates) = tips;
    colnames(allrates) = as<CharacterVector>(sample_times);

    // save a copy of these for constant-time access instead of
    // multiple pointer derefs in a hot loop
    NumericVector edge_length = as<NumericVector>(x["edge.length"]);
    NumericMatrix edge = as<NumericMatrix>(x["edge"]);
    NumericVector edge1 = edge(_, 0);
    NumericVector edge2 = edge(_, 1);
    for (std::size_t ii = 0; ii != ntips; ++ii) {
        // Although Rcpp is 0 offset, edges in the phylo edge matrix use R numbering
        std::size_t node = ii + 1;

        std::size_t count = 0;
        double cur_time = 0;
        std::size_t sample_idx = 0;

        while (node != rootnode) {
            // avoid horrible incompatible template errors by extracting a scalar value
            // instead of relying on magical Rcpp sugared subsets
            std::size_t parent = c_which(edge2, node);
            cur_time += edge_length[parent];

            while (sample_times[sample_idx] < cur_time) {
                allrates(ii, sample_idx) = count;
                sample_idx++;
            }
            node = edge1[parent];
            if (node == rootnode) {
                // add in the last count if we've reached the end of the tree but still have
                // a sampled time that hasn't been included
                while (sample_idx < sample_times.size()) {
                    allrates(ii, sample_idx) = count;
                    sample_idx++;
                }
            }
            count++;
        }

    }
    return allrates;
}
