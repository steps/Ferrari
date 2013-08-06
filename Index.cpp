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
//--------------------------------------------------------------------------------------------------
#include <assert.h>
#include <iostream>
#include <math.h>
#include <set>
#include <string.h>
#include <queue>
//--------------------------------------------------------------------------------------------------
//! Constructor.
/*!
 * \param g pointer to graph
 * \param s number of seed nodes to use for pruning
 * \param k parameter to control number of index entries
 * \param global indicator for global or local space restriction
 */
Index::Index(Graph* g, unsigned s, unsigned k, bool global) :
    g(g), n_(g->num_nodes()), s_(s), k_(k), global_(global), queryId(0), expanded(
        0) {
}
//--------------------------------------------------------------------------------------------------
//! Destructor.
Index::~Index() {
  for (unsigned i = 0; i < n_; ++i) {
    if (intervals[i]) {
      delete intervals[i];
      intervals[i] = 0;
    }
  }
}
//--------------------------------------------------------------------------------------------------
//! Build the index.
void Index::build() {
  visited = std::vector<unsigned>(n_, 0);
  expanded = 0;

  // determine seeds
  if (s_) {
    seed_nodes = std::vector<unsigned>(); // collection of seed nodes
    is_seed = std::vector<bool>(n_, false);  // indicator whether certain node is a seed 
    generate_seeds_degree(s_); // generate seed nodes by choosing nodes with highest degree

    // initialize bitmaps for seed nodes
    seed_in = std::vector<bitset>(n_, bitset(s_)); // for each node bit is set if seed reaches node
    seed_out = std::vector<bitset>(n_, bitset(s_)); // for each node bit is set if node reaches seed
    for (unsigned i = 0; i < s_; ++i) {
      // initialize for seed nodes themselves
      seed_in[seed_nodes[i]][i] = 1;
      seed_out[seed_nodes[i]][i] = 1;
    }
  }

  // filters
  tlevel_ = std::vector<unsigned>(n_, 0);  // topological level of vertex (see paper)
  torder_ = std::vector<unsigned>(n_, ~0u); // topological order number of vertex

  // postorder ids of vertices
  id_ = std::vector<unsigned>(n_, ~0u);

  // topological ordering
  ordered_nodes_ = std::vector<unsigned>();

  // parent in tree (specifies tree cover)
  parent_ = std::vector<unsigned>(n_, ~0u);

  // create index
  visit_fw();
  visit_bw();

  // intervals
  intervals = std::vector<IntervalList*>(n_, 0);

  if (global_)
    assign_sketches_global();
  else
    assign_sketches();
}
//--------------------------------------------------------------------------------------------------
//! Determine reachability using DFS traversal
/*!
 * \param x source node
 * \param y target node
 */
bool Index::reachable_dfs(unsigned x, unsigned y) {
  ++queryId; // bookkeeping for fast query processing
  return __reachable_dfs(x, y);
}
//--------------------------------------------------------------------------------------------------
//! Recursive call for DFS-based reachability processing
/*!
 * \param x source node
 * \param y target node
 */
bool Index::__reachable_dfs(unsigned x, unsigned y) {
  if (x == y)
    return true;
  if (visited[x] == queryId)
    return false; // don't visit vertex again
  visited[x] = queryId; // remember visit
  ++expanded;
  const std::vector<unsigned> *nb = g->get_neighbors(x);
  for (std::vector<unsigned>::const_iterator it = nb->begin(); it != nb->end();
      ++it) {
    // recursively query successors
    if (__reachable_dfs(*it, y)) {
      return true;
    }
  }
  return false; // x could not reach y in DFS traversal
}
//--------------------------------------------------------------------------------------------------
//! Determine reachability using BFS traversal
bool Index::reachable_bfs(unsigned x, unsigned y) {
  if (x == y)
    return true;
  ++queryId;
  std::deque<unsigned> queue(1, x); // put x in expansion queue
  unsigned v;
  const std::vector<unsigned> *nb;
  while (!queue.empty()) { // process queue while unexpanded vertices remain
    v = queue.front();
    queue.pop_front();
    if (visited[v] == queryId) 
      continue; // don't visit vertex again
    visited[v] = queryId; // mark visit
    ++expanded;
    nb = g->get_neighbors(v);
    for (std::vector<unsigned>::const_iterator it = nb->begin();
        it != nb->end(); ++it) {
      if (y == *it)
        return true;
      // add neighbor to end of queue for later expansion
      queue.push_back(*it);
    }
  }
  return false;
}
//--------------------------------------------------------------------------------------------------
//! Determin reachability from x to y using index
/*!
 * \param x source node
 * \param y target node
 */
bool Index::reachable(unsigned x, unsigned y) {
  ++queryId;
  if (!intervals[x]) {
    return x == y; // x has no interval assigned, thus no successors
  }

  if (s_) { // are we using seeds?
    if (seed_out[x].intersects(seed_in[y])) {
      // there is a seed reachable from x that can reach y, thus graph contains a path from x to y
      return true;
    }
  }
  if (x == y) {
    return true;
  }
  if (tlevel_[x] <= tlevel_[y] || torder_[x] > torder_[y]) {
    // can we use pruning techniques because x cannot reach y by topological properties?
    return false;
  }

  // identify whether y's id is outside intervals of x or inside exact/approximate intervals
  switch (intervals[x]->contains(id_[y])) {
  case IntervalList::NOT:
    // id of y is outside x's intervals (-> not reachable)
    return false;
  case IntervalList::YES:
    // id of y is inside exact interval at x (-> reachable)
    return true;
  default:
    // id of y was contained in approximate interval at x (-> potentially reachable)
    // we need to query neighbors recursively
    const std::vector<unsigned> *nb = g->get_neighbors(x);
    for (std::vector<unsigned>::const_iterator it = nb->begin();
        it != nb->end(); ++it) {
      if (y == *it
          || (visited[*it] < queryId && __reachable(*it, y))) {
        // neighbor was target or neighbor can reach target (-> reachable)
        return true;
      }
    }
    // we could not find a connection, thus unreachable
    return false;
  }
}
//--------------------------------------------------------------------------------------------------
//! Recursive call for reachability using index
/*!
 * \param x source node
 * \param y target node
 */
bool Index::__reachable(const unsigned& x, const unsigned& y) {
  visited[x] = queryId;
  ++expanded;

  // can topological properties reveal that y is not reachable from x?
  if (tlevel_[x] <= tlevel_[y] || torder_[x] > torder_[y]) {
    return false;
  }

  if (x == y) {
    return true;
  }

  // x is a leaf (-> y unreachable)
  if (!g->get_neighbors(x)) {
    return false;
  }

  // identify whether y's id is outside intervals of x or inside exact/approximate intervals
  switch (intervals[x]->contains(id_[y])) {
  case IntervalList::NOT:
    // id of y is outside x's intervals (-> not reachable)
    return false;
  case IntervalList::YES:
    // id of y is inside exact interval at x (-> reachable)
    return true;
  default:
    // id of y was contained in approximate interval at x (-> potentially reachable)
    // we need to query neighbors recursively
    const std::vector<unsigned> *nb = g->get_neighbors(x);
    for (std::vector<unsigned>::const_iterator it = nb->begin();
        it != nb->end(); ++it) {
      if (*it == y
          || (visited[*it] < queryId && __reachable(*it, y))) {
        // neighbor was target or neighbor can reach target (-> reachable)
        return true;
      }
    }
    // we could not find a connection, thus unreachable
    return false;
  }
}
//--------------------------------------------------------------------------------------------------
//! Visit vertices upwards starting from the leaves and build relevant index structures
/*!
 * Execution of this method will build the tree cover as described in the paper, as well
 * as determine the values for topological level, that are used for pruning in the query
 * processing step. Further, the bitmaps for outgoing seeds are populated.
 */
void Index::visit_bw() {
  std::deque<unsigned> leaves = *(g->get_leaves());  // leaves of the graph
  std::vector<unsigned> deg = *(g->get_degrees());   // degree of the vertices 
  const std::vector<unsigned> *pd;  // predecessors of a vertex (nodes with edge to vertex)
  unsigned v;
  while (!leaves.empty()) {
    v = leaves.front();
    leaves.pop_front();
    pd = g->get_predecessors(v);
    for (std::vector<unsigned>::const_iterator it = pd->begin();
        it != pd->end(); ++it) {
      // record topological level (see paper)
      tlevel_[*it] = std::max(tlevel_[*it], 1 + tlevel_[v]);
      if (parent_[v] == ~0u || torder_[parent_[v]] < torder_[*it]) {
        // specify parent in the tree cover by choosing parent with highest top. order (see paper)
        parent_[v] = *it;
      }

      // mark for a vertex which seeds it can reach
      if (s_)
        seed_out[*it] |= seed_out[v];
      --deg[*it];
      if (!deg[*it]) {
        // By removing this node from children of the predecessor, the
        // predecessor became a leaf? Then mark it for processing in next round
        leaves.push_back(*it);
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------
//! Visit vertices downwards starting from the roots (nodes without incoming edges)
/*!
 *
 */
void Index::visit_fw() {
  std::vector<unsigned> *roots = g->get_roots(); // the roots of the dag
  std::queue<unsigned> queue;
  for (std::vector<unsigned>::const_iterator it = roots->begin();
      it != roots->end(); ++it)
    queue.push(*it);
  std::vector<unsigned> indeg = *(g->get_indegrees());
  const std::vector<unsigned> *nb;
  unsigned v, to = 0;
  while (!queue.empty()) {
    v = queue.front();
    queue.pop();
    torder_[v] = to++; // record topological order number
    ordered_nodes_.push_back(v); // create the topological ordering
    nb = g->get_neighbors(v);
    for (std::vector<unsigned>::const_iterator it = nb->begin();
        it != nb->end(); ++it) {

      // mark for a vertex which seeds can reach it
      if (s_)
        seed_in[*it] |= seed_in[v];
      --indeg[*it];
      if (!indeg[*it]) {
        // remove vertex from child's predecessors. If it becomes a root, mark
        // it for processing in next round
        queue.push(*it);
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------
//! Assign reachability intervals to the vertices in the graph (builds actual index)
/*!
 * This is the variant FERRARI-L with a local restriction on the index size (at most k intervals
 * per node)
 */
void Index::assign_sketches() {
  // assign intervals based on tree cover specified by parent_
  // traverse dfs
  std::vector<unsigned> stack(*(g->get_roots()));
  unsigned n = g->num_nodes();
  std::vector<unsigned> mid = std::vector<unsigned>(n, ~0u);
  std::vector<unsigned> next_child(n, 0);
  unsigned _id = 0, v, c;
  std::vector<unsigned> const *nb;

  // perform a stack-based DFS traversal of the graph
  while (!stack.empty()) {
    v = stack.back();
    nb = g->get_neighbors(v);
    if (mid[v] < ~0u && next_child[v] >= nb->size()) {
      // this node can be visited, since all children have been labeled with their id
      stack.pop_back();
      id_[v] = _id++; // this is the postorder id of the node, the right end of its tree interval
      intervals[v] = new IntervalList(mid[v], id_[v]);
    } else {
      // we see node for the first time, the minimum id among the descendant must thus
      // be the current id
      mid[v] = _id;
      if (next_child[v] < nb->size()) {
        c = nb->at(next_child[v]++);
        if (parent_[c] == v) {
          // if this node is a parent of c in the tree cover, traverse along this path
          stack.push_back(c);
        }
      }
    }
  }

  // visit nodes in reverse topological order and assign the intervals
  for (std::vector<unsigned>::const_reverse_iterator it =
      ordered_nodes_.rbegin(); it != ordered_nodes_.rend(); ++it) {
    nb = g->get_neighbors(*it);
    for (std::vector<unsigned>::const_iterator nbit = nb->begin();
        nbit != nb->end(); ++nbit) {
      // this edge is a shortcut in the tree cover, we don't need to do anything,
      // since child's intervals will be propagated here along the tree path in the tree cover
      if (!(parent_[*nbit] == *it) && (mid[*it] <= id_[*nbit])
          && (id_[*nbit] < id_[*it])) {
        // forward edge
      } else {
        // a regular edge, we have to merge child's intervals into current node
        intervals[*it]->merge(*intervals[*nbit]);
      }
    }
    // ensure at most k interval by merging adjacent intervals
    intervals[*it]->restrict(k_);
  }

  // make sure that leaves don't use space for their trivial interval
  for (unsigned i = 0; i < g->num_nodes(); ++i) {
    nb = g->get_neighbors(i);
    if (!nb->size()) { // leaf node
      delete intervals[i];
      intervals[i] = 0;
    }
  }
}
//--------------------------------------------------------------------------------------------------
//!  Assign reachability intervals to the vertices in the graph (builds actual index)
/*!
 * This is the variant FERRARI-G with a global restriction on the index size (at most nk intervals
 * over all nodes)
 */
void Index::assign_sketches_global() {
  // assign intervals based on tree cover specified by parent_
  // traverse dfs
  std::vector<unsigned> stack(*(g->get_roots()));
  unsigned n = g->num_nodes();
  std::vector<unsigned> mid = std::vector<unsigned>(n, ~0u); // the minimum id among successors
  std::vector<unsigned> next_child(n, 0);
  unsigned _id = 0, v, c;
  std::vector<unsigned> const *nb;

  // perform a stack-based DFS traversal of the Graph
  while (!stack.empty()) {
    v = stack.back();
    nb = g->get_neighbors(v);
    if (mid[v] < ~0u && next_child[v] >= nb->size()) {
      // this node can be visited, since all children have been labeled with their id
      stack.pop_back();
      id_[v] = _id++; // this is the postorder id of the node, the right end of its tree interval
      intervals[v] = new IntervalList(mid[v], id_[v]);
    } else {
      // we see node for the first time, the minimum id among the descendant must thus
      // be the current id
      mid[v] = _id;
      if (next_child[v] < nb->size()) {
        c = nb->at(next_child[v]++);
        if (parent_[c] == v) {
          // if this node is a parent of c in the tree cover, traverse along this path
          stack.push_back(c);
        } 
      }
    }
  }

  // parameters for global index size
  unsigned multiplier = 4;  // on first shot, restrict interval list to 4 times the average budget
  unsigned leaf_count = g->get_leaves()->size();
  unsigned budget = k_ * (n + leaf_count) / n; // average budget a node has

  /* maximum space we can allocate, assuming leaves don't consume space for intervals */
  unsigned max_space = n * k_ + leaf_count; 
  unsigned current_space = 0; // number of currently assigned intervals
  std::vector<unsigned> restriction_queue; // queue of nodes that can be restricted
  restriction_queue.reserve(n);
  std::vector<unsigned>* deg = g->get_degrees();
  degreecompare comp(deg); // comparator ensuring nodes with lower degree get restricted first

  // visit nodes in reverse topological order and assign intervals
  for (std::vector<unsigned>::const_reverse_iterator it =
      ordered_nodes_.rbegin(); it != ordered_nodes_.rend(); ++it) {
    nb = g->get_neighbors(*it);
    for (std::vector<unsigned>::const_iterator nbit = nb->begin();
        nbit != nb->end(); ++nbit) {
      // this edge is a shortcut in the tree cover, we don't need to do anything,
      // since child's intervals will be propagated here along the tree path in the tree cover
      if (!(parent_[*nbit] == *it) && (mid[*it] <= id_[*nbit])
          && (id_[*nbit] < id_[*it])) {
        // forward edge
      } else {
        // a regular edge, we have to merge child's intervals into current nodes
        intervals[*it]->merge(*intervals[*nbit]);
      }
    }
    // first restrict interval list with a generous budget
    intervals[*it]->restrict(multiplier * budget);
    current_space += intervals[*it]->size(); // update space consumption
    unsigned v;

    // are we over budget? then start restricting nodes in ascending order of degree
    while (current_space > max_space) {
      if (!restriction_queue.empty()) {
        std::pop_heap(restriction_queue.begin(), restriction_queue.end(), comp);
        v = restriction_queue.back();
        restriction_queue.pop_back();

        // we have freed some space
        current_space -= (intervals[v]->size() - budget);

        // restrict node to its strict budget
        intervals[v]->restrict(budget);
      } else {
        // we could not restrict more nodes with lower degree, so current node must be restricted
        if (intervals[*it]->size() > budget) {
          current_space -= intervals[*it]->size() - budget;
          intervals[*it]->restrict(budget);
        }
      }
    }

    // record node for possible restriction in the future
    if (intervals[*it]->size() > budget) {
      restriction_queue.push_back(*it);
      std::push_heap(restriction_queue.begin(), restriction_queue.end(), comp);
    }
  }

  // free space used by the trivial leaf intervals
  for (unsigned i = 0; i < g->num_nodes(); ++i) {
    nb = g->get_neighbors(i);
    if (!nb->size()) { // leaf node
      delete intervals[i];
      intervals[i] = 0;
    }
  }
}
//--------------------------------------------------------------------------------------------------
//! Generate random vertices as seeds
/*!
 * \param s number of seed vertices to use
 */
void Index::generate_seeds_random(unsigned s) {
  std::srand(std::time(0));
  std::set<unsigned> seed_set;
  unsigned i, n = g->num_nodes();
  while (seed_set.size() < s) {
    i = std::rand() % n + 1;
    // to qualify for a seed, node must have an incoming and outgoing edge
    if (g->get_indegrees()->at(i) && g->get_degrees()->at(i)) {
      seed_set.insert(i);
      is_seed[i] = true;
    }
  }
  std::copy(seed_set.begin(), seed_set.end(), std::back_inserter(seed_nodes));
}
//--------------------------------------------------------------------------------------------------
//! Pick nodes with highest degrees as seeds
/*!
 * \param s number of seed vertices to use
 */
void Index::generate_seeds_degree(unsigned s) {
  // select nodes with maximum degree as seeds
  std::vector<unsigned>* deg = g->get_degrees();
  //std::vector<unsigned>* indeg = g->get_indegrees();
  degreecompare comp(deg);
  //combined_degreecompare comp(deg, indeg);
  for (unsigned i = 0; i < g->num_nodes(); ++i) {
    if (seed_nodes.size() < s && deg->at(i) > 1) {
      seed_nodes.push_back(i);
      std::push_heap(seed_nodes.begin(), seed_nodes.end(), comp);
    } else {
      if (/*deg->at(seed_nodes.front()) * */deg->at(i) > 1
          && deg->at(seed_nodes.front()) < /*deg->at(i) * */deg->at(i)) {
        std::pop_heap(seed_nodes.begin(), seed_nodes.end(), comp);
        seed_nodes.pop_back();
        seed_nodes.push_back(i);
        std::push_heap(seed_nodes.begin(), seed_nodes.end(), comp);
      }
    }
  }
  for (std::vector<unsigned>::const_iterator it = seed_nodes.begin();
      it != seed_nodes.end(); ++it) {
    is_seed[*it] = true;
  }
  s_ = std::min((unsigned) seed_nodes.size(), s_);
}
//--------------------------------------------------------------------------------------------------
//! Reset bookkeeping data structures used for query processing 
unsigned Index::reset() {
  unsigned _expanded = expanded;
  expanded = 0;
  queryId = 0;
  memset(&visited[0], 0, sizeof(visited[0]) * visited.size());
  return _expanded;
}
//--------------------------------------------------------------------------------------------------
//! Compute the path from x to y (for debugging purposes)
/*!
 * \param x source node
 * \param y target node
 * \param p pointer to vector that will hold the resulting path
 * \return true if path exists, false if not
 */
bool Index::path(unsigned x, unsigned y, std::vector<unsigned>* p) {
  if (x == y) {
    p->push_back(x);
    return true;
  } else {
    if (!reachable_dfs(x, y)) {
      return false;
    } else {
      p->push_back(x);
      const std::vector<unsigned> *nb = g->get_neighbors(x);
      for (std::vector<unsigned>::const_iterator it = nb->begin();
          it != nb->end(); ++it) {
        if (reachable_dfs(*it, y)) {
          return path(*it, y, p);
        }
      }
    }
  }
  return true;
}
//--------------------------------------------------------------------------------------------------
