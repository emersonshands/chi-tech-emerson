#include "chi_directed_graph.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include <sstream>
#include <algorithm>

//###################################################################
/** Adds a vertex to the graph with a supplied id.*/
void chi_graph::DirectedGraph::
  VertexAccessor::AddVertex(size_t id, void* context)
{
  vertices.emplace_back(id, context);
  vertex_valid_flags.push_back(true);
}

/** Adds a vertex to the graph where the ID is assigned to
 * the number of vertices already loaded on the graph.
 * For example, if there are 3 vertices on the graph (with
 * IDs 0 through 2) then the next vertex (this one) will
 * be assigned and ID of 3.*/
void chi_graph::DirectedGraph::
VertexAccessor::AddVertex(void* context)
{
  vertices.emplace_back(vertices.size(), context);
  vertex_valid_flags.push_back(true);
}


//###################################################################
/** Removes a vertex from the graph.*/
void chi_graph::DirectedGraph::
  VertexAccessor::RemoveVertex(size_t v)
{
  ChiLogicalError(v>=vertices.size(),"Error removing vertex.")

  auto& vertex = vertices[v];

  //=================================== Get adjacent vertices
  auto num_us = vertex.us_edge.size();
  auto num_ds = vertex.ds_edge.size();
  std::vector<size_t> adj_verts;
  adj_verts.reserve(num_us+num_ds);

  for (size_t u : vertex.us_edge)
    adj_verts.push_back(u);

  for (size_t u : vertex.ds_edge)
    adj_verts.push_back(u);

  //=================================== Remove v from all u
  for (size_t u : adj_verts)
  {
    vertices[u].us_edge.erase(v);
    vertices[u].ds_edge.erase(v);
  }

  vertex_valid_flags[v] = false;
}

//###################################################################
/** Accesses a vertex from the graph.*/
chi_graph::GraphVertex& chi_graph::DirectedGraph::
  VertexAccessor::operator[](size_t v)
{
  if (not vertex_valid_flags[v])
    chi::log.LogAllError()
      << "chi_graph::DirectedGraph::VertexAccessor: "
         "Invalid vertex accessed. Vertex may have been removed.";
  return vertices[v];
}

//###################################################################
/** Adds a vertex to the graph. By default <I>context</I> is
 * assumed to be nullptr.*/
void chi_graph::DirectedGraph::AddVertex(size_t id,void* context/*=nullptr*/)
{
  vertices.AddVertex(id, context);
}

//###################################################################
/** Adds a vertex to the graph. By default <I>context</I> is
 * assumed to be nullptr and <I>id</I> is assumed to be assigned
 * automatically. In
 * the latter case the vertex id will be the same as the order
 * in which it was added (0,1,2,3,etc ... will have id's 0,1,2,3,etc)*/
void chi_graph::DirectedGraph::AddVertex(void* context/*=nullptr*/)
{
  vertices.AddVertex(context);
}

//###################################################################
/** Removes a vertex from the graph. This method does not
 * free any context related data.*/
void chi_graph::DirectedGraph::
  RemoveVertex(size_t v)
{
  vertices.RemoveVertex(v);
}

//###################################################################
/** Adds an edge to the graph. Range checks are supplied by the
 * vertex accessor.*/
bool chi_graph::DirectedGraph::
  AddEdge(size_t from, size_t to, double weight)
{
  vertices[from].ds_edge.insert(to);
  vertices[to].us_edge.insert(from);

  vertices[from].ds_weights[to] = weight;
  vertices[to].us_weights[from] = weight;

  return true;
}

//###################################################################
/**Remove an edge from the graph. Range checks are supplied by the
 * vertex accessor.*/
void chi_graph::DirectedGraph::RemoveEdge(size_t from, size_t to)
{
  vertices[from].ds_edge.erase(to);
  vertices[to].us_edge.erase(from);
}

//###################################################################
/** Depth-First-Search main recursive algorithm. This is the recursive
 * portion of the method below this one
 * (chi_graph::DirectedGraph::DepthFirstSearch).*/
void chi_graph::DirectedGraph::
  DFSAlgorithm(std::vector<size_t> &traversal,
               std::vector<bool> &visited,
               size_t cur_vid)
{
  traversal.push_back(cur_vid);
  visited[cur_vid] = true;

  for (auto v : vertices[cur_vid].ds_edge)
    if (not visited[v])
      DFSAlgorithm(traversal,visited,v);

}

//###################################################################
/**SCC main recursive algorithm. This is the recursive call for the
 * method defined below this one
 * (chi_graph::DirectedGraph::FindStronglyConnectedConnections).*/
void chi_graph::DirectedGraph::SCCAlgorithm(
  size_t u, int& time,
  std::vector<int>& disc,
  std::vector<int>& low,
  std::vector<bool>& on_stack,
  std::stack<size_t>& stack,
  std::vector<std::vector<size_t>>& SCCs)
{
  // Init discovery time and low value
  disc[u] = low[u] = ++time;
  stack.push(u);
  on_stack[u] = true;

  for (auto v : vertices[u].ds_edge)
  {
    if (disc[v] == -1)
    {
      SCCAlgorithm(v,time,disc,low,on_stack,stack,SCCs);
      low[u] = std::min(low[u],low[v]);
    }
    else if (on_stack[v])
      low[u] = std::min(low[u],disc[v]);
  }

  size_t w=0;
  if (low[u] == disc[u])
  {
    std::vector<size_t> sub_SCC;
    while (stack.top() != u)
    {
      w = stack.top();
      sub_SCC.push_back(w);
      on_stack[w] = false;
      stack.pop();
    }
    w = stack.top();
    sub_SCC.push_back(w);
    if (sub_SCC.size() > 1) SCCs.push_back(sub_SCC);
    on_stack[w] = false;
    stack.pop();
  }

}

//###################################################################
/**Find strongly connected components. This method is the implementation
 * of Tarjan's algorithm [1].
 *
 * [1] Tarjan R.E. "Depth-first search and linear graph algorithms",
 *     SIAM Journal on Computing, 1972.
 *
 * It returns collections of vertices that form strongly connected
 * components excluding singletons.*/
std::vector<std::vector<size_t>> chi_graph::DirectedGraph::
  FindStronglyConnectedComponents()
{
  size_t V = vertices.size();

  std::vector<int>     disc(V,-1);        // Discovery times
  std::vector<int>     low(V,-1);         // Earliest visited vertex
  std::vector<bool>    on_stack(V,false); // On stack flags
  std::stack<size_t> stack;             // Stack

  std::vector<std::vector<size_t>> SCCs;  // Collection of SCCs

  int time = 0;

  for (size_t v=0; v<V; ++v)
    if (disc[v]==-1)
      SCCAlgorithm(v,time,disc,low,on_stack,stack,SCCs);


  return SCCs;
}


//###################################################################
/** Generates a topological sort. This method is the implementation
 * of Kahn's algorithm [1].
 *
 * [1] Kahn, Arthur B. (1962), "Topological sorting of large networks",
 *     Communications of the ACM, 5 (11): 558–562
 *
 * \return Returns the vertex ids sorted topologically. If this
 *         vector is empty the algorithm failed because it detected
 *         cyclic dependencies.*/
std::vector<size_t> chi_graph::DirectedGraph::GenerateTopologicalSort()
{
  bool has_cycles=false;
  std::vector<size_t> L;
  std::vector<GraphVertex*> S;

  L.reserve(vertices.size());
  S.reserve(vertices.size());

  //======================================== Make a copy of the graph
  auto cur_vertices = vertices;

  //======================================== Identify vertices that
  //                                         have no incoming edge
  for (auto& vertex : cur_vertices)
    if (vertex.us_edge.empty())
      S.push_back(&vertex);

  if (S.empty())
  {
    has_cycles = true;
    goto endofalgo;
  }

  //======================================== Repeatedly remove
  //                                         vertices
  while (not S.empty())
  {
    GraphVertex* node_n = S.back(); size_t n=node_n->id;
    S.erase(S.end()-1);

    L.push_back(n);
    auto nodes_m = node_n->ds_edge;
    for (size_t m : nodes_m)
    {
      GraphVertex* node_m = &cur_vertices[m];

      //remove edge from graph
      node_n->ds_edge.erase(m);
      node_m->us_edge.erase(n);

      if (node_m->us_edge.empty())
        S.push_back(node_m);
    }
  }

  endofalgo:
  {
    if (has_cycles)
      return {};
    if (L.size() != vertices.size())
      return {};
  }

  return L;
}

//###################################################################
/**Finds a sequence that minimizes the Feedback Arc Set (FAS). This
 * algorithm implements the algorithm depicted in [1].
 *
 * [1] Eades P., Lin X., Smyth W.F., "Fast & Effective heuristic for
 *     the feedback arc set problem", Information Processing Letters,
 *     Volume 47. 1993.*/
std::vector<size_t> chi_graph::DirectedGraph::
  FindApproxMinimumFAS()
{
  auto GetVertexDelta = [](chi_graph::GraphVertex& vertex)
  {
    double delta = 0.0;
    for (size_t ds : vertex.ds_edge)
      delta += 1.0*vertex.ds_weights[ds];

    for (size_t us : vertex.us_edge)
      delta -= 1.0*vertex.us_weights[us];

    return delta;
  };

  auto& TG = *this;

  //==================================== Execute GR-algorithm
  std::vector<size_t> s1,s2,s;
  while (TG.vertices.GetNumValid()>0)
  {
    //======================== Remove sinks
    while (TG.GetNumSinks()>0)
    {
      for (auto& u : TG.vertices)
        if (u.ds_edge.empty())
        {
          TG.RemoveVertex(u.id);
          s2.push_back(u.id);
          break;
        }
    }//G contains sinks

    //======================== Remove sources
    while (TG.GetNumSources()>0)
    {
      for (auto& u : TG.vertices)
        if (u.us_edge.empty())
        {
          TG.RemoveVertex(u.id);
          s1.push_back(u.id);
          break;
        }
    }//G contains sinks

    //======================== Get max delta
    std::pair<int,double> max_delta(-1,-100.0);
    for (auto& u : TG.vertices)
    {
      double delta = GetVertexDelta(u);
      if (delta > max_delta.second)
        max_delta = std::make_pair(u.id,delta);
    }

    //======================== Remove max delta
    TG.RemoveVertex(max_delta.first);
    s1.push_back(max_delta.first);
  }


  //========================== Make appr. minimum FAS sequence
  s.reserve(s1.size() + s2.size());
  for (size_t u : s1) s.push_back(u);
  for (size_t u : s2) s.push_back(u);

  return s;
}


//###################################################################
/**Prints the graph in Graphviz format.*/
void chi_graph::DirectedGraph::PrintGraphviz(int location_mask)
{
  if (chi::mpi.location_id != location_mask) return;

  std::stringstream o;
  std::string offset("    ");
  o << "Printing directed graph:\n";
  o << "digraph DG {\n";

  o << offset << "splines=\"FALSE\";\n";
  o << offset << "rankdir=\"LR\";\n\n";


  o << offset << "/* Vertices */\n";
  for (auto& v : vertices)
    o << offset << v.id << " [shape=\"circle\"]\n";

  o << "\n" << offset << "/* Edges */\n";
  for (auto& v : vertices)
    for (int w : v.ds_edge)
      o << offset << v.id << " -> " << w << "\n";

  o << "}\n";

  std::cout << o.str();
}

//###################################################################
/**Prints a sub-graph in Graphviz format.*/
void chi_graph::DirectedGraph::
  PrintSubGraphviz(const std::vector<int>& verts_to_print,
                   int location_mask)
{
  if (chi::mpi.location_id != location_mask) return;

  std::stringstream o;
  std::string offset("    ");
  o << "Printing directed graph:\n";
  o << "digraph DG {\n";

  o << offset << "splines=\"FALSE\";\n";
  o << offset << "rankdir=\"LR\";\n\n";


  o << offset << "/* Vertices */\n";
  for (int v : verts_to_print)
    o << offset << v << " [shape=\"circle\"]\n";

  o << "\n" << offset << "/* Edges */\n";
  for (int v : verts_to_print)
    for (int w : vertices[v].ds_edge)
    {
      if (std::find(verts_to_print.begin(),
                    verts_to_print.end(),
                    w) != verts_to_print.end())
        o << offset << v << " -> " << w << "\n";
    }

  o << "}\n";

  std::cout << o.str();
}

//###################################################################
std::vector<std::pair<size_t,size_t>>
chi_graph::DirectedGraph::RemoveCyclicDependencies()
{
  std::vector<std::pair<size_t,size_t>> edges_to_remove;

  //============================================= Utility lambdas
  auto IsInList = [](std::vector<size_t>& list, size_t val)
  {
    return std::find(list.begin(),list.end(),val) != list.end();
  };


  //============================================= Find initial SCCs
  auto SCCs = FindStronglyConnectedComponents();

  int iter=0;
  while (not SCCs.empty())
  {
    if (chi::log.GetVerbosity() >= chi_objects::ChiLog::LOG_LVL::LOG_0VERBOSE_2)
      chi::log.LogAll()
        << "Inter cell cyclic dependency removal. Iteration " << ++iter;

    //============================================= Remove bi-connected then
    //                                              tri-connected SCCs then
    //                                              n-connected
    for (auto& subDG : SCCs)
    {
      //====================================== If bi-connected
      if (subDG.size()==2)
      {
        RemoveEdge(subDG.front(), subDG.back());
        edges_to_remove.emplace_back(subDG.front(), subDG.back());
      }//bi-connected
        //====================================== If tri-connected
      else if (subDG.size()==3)
      {
        bool found=false;
        for (size_t u : subDG)
        {
          for (size_t v : vertices[u].ds_edge)
            if (IsInList(subDG,v))
            {
              found=true;
              RemoveEdge(u, v);
              edges_to_remove.emplace_back(u, v);
              break;
            }
          if (found) break;
        }//for u
      }//tri-connected
        //====================================== If n-connected
      else
      {
        //=============================== Add vertices to temporary graph
        chi_graph::DirectedGraph TG; //Temp Graph
        for (size_t k=0; k<subDG.size(); ++k)
          TG.AddVertex();

        //=============================== Add local connectivity
        int mapping_u = 0;
        for (auto u : subDG)
        {
          for (auto v : vertices[u].ds_edge)
          {
            auto mapv = std::find(subDG.begin(),subDG.end(),v);
            if (mapv != subDG.end())
            {
              size_t mapping_v = mapv - subDG.begin();
              TG.AddEdge(mapping_u,mapping_v,vertices[u].ds_weights[v]);
            }
          }//for v

          ++mapping_u;
        }//for u

        //=============================== Make a copy of the graph verts
        std::vector<chi_graph::GraphVertex> verts_copy;
        verts_copy.reserve(TG.vertices.size());
        for (auto& v : TG.vertices)
          verts_copy.push_back(v);

        //=============================== Solve the minimum Feedback
        //                                     Arc Set (FAS) problem
        auto s = TG.FindApproxMinimumFAS();

        //========================== Build a sequence map
        // This maps original sequence
        // to minFAS sequence. i.e. originally
        // we had v=0,1,2,3... and afterwards we
        // something like s=7,3,1,5,0,....
        // smap[v] then gives the position of v in s
        std::vector<int> smap(s.size(),-1);
        int count=0;
        for (size_t u: s)
          smap[u] = count++;

        //========================== Build edges to remove
        std::vector<std::pair<int,int>> edges_to_rem;
        for (auto& u : verts_copy)
        {
          int cur_map = smap[u.id];
          for (size_t v : u.ds_edge)
          {
            int adj_map = smap[v];
            if (adj_map < cur_map)
              edges_to_rem.emplace_back(u.id,v);
          }
        }

        for (auto& edge : edges_to_rem)
        {
          size_t u = subDG[edge.first];
          size_t v = subDG[edge.second];
          RemoveEdge(u, v);
          edges_to_remove.emplace_back(u, v);
        }

      }//n-connected
    }//for sub-DG


    //============================================= Find SSCs again
    // This step is like an insurance policy for if
    // something came through. There should be no SSCs
    // after the minFAS process, however, we just look
    // again in-case.
    SCCs = FindStronglyConnectedComponents();
  }

  return edges_to_remove;
}

//###################################################################
/**Clears all the data structures associated with the graph.*/
void chi_graph::DirectedGraph::Clear()
{
  vertices.clear();
}


//###################################################################
/**Destructor.*/
chi_graph::DirectedGraph::~DirectedGraph()
{
  Clear();
}
