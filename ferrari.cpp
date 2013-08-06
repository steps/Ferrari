//--------------------------------------------------------------------------------------------------
// Ferrari Reachability Index
// (c) 2012 Stephan Seufert. Web site: http://www.mpi-inf.mpg.de/~sseufert
//
// This work is licensed under the Creative Commons
// Attribution-Noncommercial-Share Alike 3.0 Unported License. To view a copy
// of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/
// or send a letter to Creative Commons, 171 Second Street, Suite 300,
// San Francisco, California, 94105, USA.
//--------------------------------------------------------------------------------------------------
#include "Index.h"
#include "Graph.h"
#include "IntervalList.h"
#include "Timer.h"
//--------------------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
//--------------------------------------------------------------------------------------------------
void print_usage() {
  std::cout << "Ferrari Reachability Index" << std::endl;
  std::cout << "Usage: ferrari -g <GRAPHFILE> -q <QUERYFILE> -k <SIZE> [-s <SEEDS> ] [-L]" 
            << std::endl;
}
//--------------------------------------------------------------------------------------------------
void read_queries(const std::string& query_file,
    std::vector<std::pair<unsigned, unsigned> > *queries) {
  std::ifstream qf(query_file.c_str(), std::ios::in);
  if (qf.is_open()) {
    std::string line;
    unsigned s, t;
    if (!qf.eof()) {
      while (qf.good()) {
        getline(qf, line);
        if (line.length()) {
          std::istringstream iss(line, std::istringstream::in);
          iss >> s;
          iss >> t;
          queries->push_back(std::make_pair(s, t));
        }
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  std::string graph_file = "", query_file = "";
  unsigned seeds = 0, k = ~0u; 
  bool global = true;
  for (int i = 0; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-g") {         // graph file
      graph_file = argv[++i];
    } else if (arg == "-s") {  // number of seeds
      seeds = atoi(argv[++i]);
    } else if (arg == "-q") {  // query file
      query_file = argv[++i];
    } else if (arg == "-k") {  // index size constraint
      k = atoi(argv[++i]);
    } else if (arg == "-L") {  // local size buget
      global = false;
    } else if (arg == "-G") {  // global size budget
      global = true;
    }
  }

  if (graph_file == "" || query_file == "" || k==~0u) {
    print_usage();
    return 1;
  }  

  std::cout << "graph file: " << graph_file << std::endl;
  std::cout << "query file: " << query_file << std::endl;
  std::cout << "number of seeds: " << seeds << std::endl;
  if (k < ~0u) {
    std::cout << "size constraint: " << k << std::endl;
  } else {
    std::cout << "size constraint: none" << std::endl;
  }

  // parse graph
  Graph *g = new Graph(graph_file);

  // build index
  Timer t1; t1.start();
  Index bm(g, seeds, k, global);
  bm.build();
  double t_index = t1.stop();
  std::cout << "assigned bitmaps (" << t_index << " ms)" << std::endl;

  // assess index size
  unsigned count = 0;
  for (unsigned i = 0; i < g->num_nodes(); ++i) {
     if (bm.get_intervals(i)) {
       count += bm.get_intervals(i)->size();
     }
  }
  std::cout << "assigned " << count << " intervals" << std::endl;
  
  // extract queries
  std::vector<std::pair<unsigned, unsigned> > queries;
  read_queries(query_file, &queries);
  std::cout << "running queries" << std::endl;
  unsigned reachable = 0;

  // probe reachability index
  t1.start();
  reachable = 0;
  for (std::vector<std::pair<unsigned, unsigned> >::const_iterator it =
      queries.begin(); it != queries.end(); ++it) {
    if (bm.reachable(it->first, it->second)) {  // reachability query
      ++reachable;
    }
  }
  double t_query = t1.stop();
  std::cout << "query processing time (" << t_query << " ms)" << std::endl;
  std::cout << bm.reset() << " expanded nodes" << std::endl;
  std::cout << reachable << "/" << queries.size() << " reachable" << std::endl;

  // calculate index size
  unsigned interval_space = count * 4 * 2 + count;  // intervals + exactness flag
  unsigned seed_space = (bm.used_seed_count() / 8) * 2 * g->num_nodes();
  unsigned idspace = g->num_nodes() * 4;
  unsigned filter_space = g->num_nodes() * 4 + g->num_nodes() * 4; // top order + top level
  std::cout << "Index Size: " << interval_space + seed_space + idspace + filter_space
            << " bytes" << std::endl;

  /* baseline methods: BFS & DFS */
  /*
  // dfs
  t1.start();
  reachable = 0;
  for (std::vector<std::pair<unsigned, unsigned> >::const_iterator it =
      queries.begin(); it != queries.end(); ++it) {
    if (bm.reachable_dfs(it->first, it->second)) {
      ++reachable;
    }
  }
  double t_dfs = t1.stop();
  std::cout << "query processing time (dfs) (" << t_dfs << " ms)" << std::endl;
  std::cout << bm.reset() << " expanded nodes" << std::endl;
  std::cout << reachable << "/" << queries.size() << " reachable" << std::endl;

  // bfs
  t1.start();
  reachable = 0;
  for (std::vector<std::pair<unsigned, unsigned> >::const_iterator it =
      queries.begin(); it != queries.end(); ++it) {
    if (bm.reachable_bfs(it->first, it->second)) {
      ++reachable;
    }
  }
  double t_bfs = t1.stop();
  std::cout << "query processing time (bfs) (" << t_bfs << " ms)" << std::endl;
  std::cout << bm.reset() << " expanded nodes" << std::endl;
  std::cout << reachable << "/" << queries.size() << " reachable" << std::endl;
  */

  delete g;
}
//--------------------------------------------------------------------------------------------------
