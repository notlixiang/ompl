/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2013, Willow Garage
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Author: Ioan Sucan, Ryan Luna */

#include "ompl/geometric/planners/prm/LazyPRM.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/geometric/planners/prm/ConnectionStrategy.h"
#include "ompl/tools/config/SelfConfig.h"
#include <boost/graph/astar_search.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/foreach.hpp>
#include <queue>

#include <iostream>
#include <fstream>

#include "GoalVisitor.hpp"

#define foreach BOOST_FOREACH

namespace ompl
{
    namespace magic
    {
        /** \brief The number of nearest neighbors to consider by
            default in the construction of the PRM roadmap */
        static const unsigned int DEFAULT_NEAREST_NEIGHBORS_LAZY = 5;

        /** \brief When optimizing solutions with lazy planners, this is the minimum
            number of path segments to add before attempting a new optimized solution
            extraction */
        static const unsigned int MIN_ADDED_SEGMENTS_FOR_LAZY_OPTIMIZATION = 5;
    }  // namespace magic
}  // namespace ompl

ompl::geometric::LazyPRM::LazyPRM(const base::SpaceInformationPtr &si, bool starStrategy)
  : base::Planner(si, "LazyPRM")
  , starStrategy_(starStrategy)
  , userSetConnectionStrategy_(false)
  , maxDistance_(0.0)
  , indexProperty_(boost::get(boost::vertex_index_t(), g_))
  , stateProperty_(boost::get(vertex_state_t(), g_))
  , weightProperty_(boost::get(boost::edge_weight, g_))
  , vertexComponentProperty_(boost::get(vertex_component_t(), g_))
  , vertexValidityProperty_(boost::get(vertex_flags_t(), g_))
  , edgeValidityProperty_(boost::get(edge_flags_t(), g_))
  , componentCount_(0)
  , bestCost_(std::numeric_limits<double>::quiet_NaN())
  , iterations_(0)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.approximateSolutions = false;
    specs_.optimizingPaths = true;

    Planner::declareParam<double>("range", this, &LazyPRM::setRange, &LazyPRM::getRange, "0.:1.:10000.");
    if (!starStrategy_)
        Planner::declareParam<unsigned int>("max_nearest_neighbors", this, &LazyPRM::setMaxNearestNeighbors,
                                            std::string("8:1000"));

    addPlannerProgressProperty("iterations INTEGER", std::bind(&LazyPRM::getIterationCount, this));
    addPlannerProgressProperty("best cost REAL", std::bind(&LazyPRM::getBestCost, this));
    addPlannerProgressProperty("milestone count INTEGER", std::bind(&LazyPRM::getMilestoneCountString, this));
    addPlannerProgressProperty("edge count INTEGER", std::bind(&LazyPRM::getEdgeCountString, this));
}

ompl::geometric::LazyPRM::~LazyPRM()
{
}

void ompl::geometric::LazyPRM::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);

    if (!nn_)
    {
        nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
        nn_->setDistanceFunction(
            std::bind(&LazyPRM::distanceFunction, this, std::placeholders::_1, std::placeholders::_2));
    }
    if (!connectionStrategy_)
    {
        if (starStrategy_)
            connectionStrategy_ =
                KStarStrategy<Vertex>(std::bind(&LazyPRM::milestoneCount, this), nn_, si_->getStateDimension());
        else
            connectionStrategy_ = KBoundedStrategy<Vertex>(magic::DEFAULT_NEAREST_NEIGHBORS_LAZY, maxDistance_, nn_);
    }
    if (!connectionFilter_)
        connectionFilter_ = [](const Vertex &, const Vertex &) { return true; };

    // Setup optimization objective
    //
    // If no optimization objective was specified, then default to
    // optimizing path length as computed by the distance() function
    // in the state space.
    if (pdef_)
    {
        if (pdef_->hasOptimizationObjective())
            opt_ = pdef_->getOptimizationObjective();
        else
        {
            opt_.reset(new base::PathLengthOptimizationObjective(si_));
            if (!starStrategy_)
                opt_->setCostThreshold(opt_->infiniteCost());
        }
    }
    else
    {
        OMPL_INFORM("%s: problem definition is not set, deferring setup completion...", getName().c_str());
        setup_ = false;
    }

    sampler_ = si_->allocStateSampler();
}

void ompl::geometric::LazyPRM::setRange(double distance)
{
    maxDistance_ = distance;
    if (!userSetConnectionStrategy_)
        connectionStrategy_ = ConnectionStrategy();
    if (isSetup())
        setup();
}

void ompl::geometric::LazyPRM::setMaxNearestNeighbors(unsigned int k)
{
    if (starStrategy_)
        throw Exception("Cannot set the maximum nearest neighbors for " + getName());
    if (!nn_)
    {
        nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
        nn_->setDistanceFunction(
            std::bind(&LazyPRM::distanceFunction, this, std::placeholders::_1, std::placeholders::_2));
    }
    if (!userSetConnectionStrategy_)
        connectionStrategy_ = ConnectionStrategy();
    if (isSetup())
        setup();
}

void ompl::geometric::LazyPRM::setProblemDefinition(const base::ProblemDefinitionPtr &pdef)
{
    Planner::setProblemDefinition(pdef);
    clearQuery();
}

void ompl::geometric::LazyPRM::clearQuery()
{
    startM_.clear();
    goalM_.clear();
    pis_.restart();
}

void ompl::geometric::LazyPRM::clear()
{
    Planner::clear();
    freeMemory();
    if (nn_)
        nn_->clear();
    clearQuery();

    componentCount_ = 0;
    iterations_ = 0;
    bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());

    addGeneratdMilestones();
}

void ompl::geometric::LazyPRM::freeMemory()
{
    foreach (Vertex v, boost::vertices(g_))
        si_->freeState(stateProperty_[v]);
    g_.clear();
}

ompl::geometric::LazyPRM::Vertex ompl::geometric::LazyPRM::addMilestone(base::State *state)
{
    Vertex m = boost::add_vertex(g_);
    stateProperty_[m] = state;
    vertexValidityProperty_[m] = VALIDITY_UNKNOWN;
    unsigned long int newComponent = componentCount_++;
    vertexComponentProperty_[m] = newComponent;
    componentSize_[newComponent] = 1;

    // Which milestones will we attempt to connect to?

    //

    const std::vector<Vertex> &neighbors = connectionStrategy_(m);  // 找到邻近的顶点
    foreach (Vertex n, neighbors)
        if (connectionFilter_(m, n))
        {
            const base::Cost weight = opt_->motionCost(stateProperty_[m], stateProperty_[n]);
            const Graph::edge_property_type properties(weight);
            const Edge &e = boost::add_edge(m, n, properties, g_).first;
            edgeValidityProperty_[e] = VALIDITY_UNKNOWN;
            uniteComponents(m, n);
        }

    nn_->add(m);

    return m;
}

ompl::base::PlannerStatus ompl::geometric::LazyPRM::solve(const base::PlannerTerminationCondition &ptc)
{
    clock_t begin_all = clock();
    checkValidity();

    //初始化初始与目标姿态
    base::GoalSampleableRegion *goal = dynamic_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());

    if (!goal)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        saveLogToFile(-1, -1);
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }

    // Add the valid start states as milestones
    while (const base::State *st = pis_.nextStart())
        startM_.push_back(addMilestone(si_->cloneState(st)));

    if (startM_.size() == 0)
    {
        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
        saveLogToFile(-1, -1);
        return base::PlannerStatus::INVALID_START;
    }

    if (!goal->couldSample())
    {
        OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
        saveLogToFile(-1, -1);
        return base::PlannerStatus::INVALID_GOAL;
    }

    // Ensure there is at least one valid goal state
    if (goal->maxSampleCount() > goalM_.size() || goalM_.empty())
    {
        const base::State *st = goalM_.empty() ? pis_.nextGoal(ptc) : pis_.nextGoal();
        if (st)
            goalM_.push_back(addMilestone(si_->cloneState(st)));

        if (goalM_.empty())
        {
            OMPL_ERROR("%s: Unable to find any valid goal states", getName().c_str());
            saveLogToFile(-1, -1);
            return base::PlannerStatus::INVALID_GOAL;
        }
    }

    unsigned long int nrStartStates = boost::num_vertices(g_);
    OMPL_INFORM("%s: Starting planning with %lu states already in datastructure", getName().c_str(), nrStartStates);

    bestCost_ = opt_->infiniteCost();
    base::State *workState = si_->allocState();
    std::pair<std::size_t, std::size_t> startGoalPair;
    base::PathPtr bestSolution;
    bool fullyOptimized = false;
    bool someSolutionFound = false;
    unsigned int optimizingComponentSegments = 0;

    // Grow roadmap in lazy fashion -- add vertices and edges without checking validity
    while (ptc == false)
    {
        ++iterations_;
        sampler_->sampleUniform(workState);
        Vertex addedVertex = addMilestone(si_->cloneState(workState));
        //找是否有已连接
        const long int solComponent = solutionComponent(&startGoalPair);
        // If the start & goal are connected and we either did not find any solution
        // so far or the one we found still needs optimizing and we just added an edge
        // to the connected component that is used for the solution, we attempt to
        // construct a new solution.

        //有未确认的解且(求解未结束或新加的点在上述解的分块上)
        if (solComponent != -1 &&
            (!someSolutionFound || (long int)vertexComponentProperty_[addedVertex] == solComponent))
        {
            // If we already have a solution, we are optimizing. We check that we added at least
            // a few segments to the connected component that includes the previously found
            // solution before attempting to construct a new solution.

            if (someSolutionFound)
            {
                if (++optimizingComponentSegments < magic::MIN_ADDED_SEGMENTS_FOR_LAZY_OPTIMIZATION)
                    continue;
                optimizingComponentSegments = 0;
            }
            Vertex startV = startM_[startGoalPair.first];
            Vertex goalV = goalM_[startGoalPair.second];
            base::PathPtr solution;
            do
            {
                solution = constructSolution(startV, goalV);
            } while (!solution && vertexComponentProperty_[startV] == vertexComponentProperty_[goalV]);
            //尝试求解直到初末姿态不在同一分块上
            if (solution)
            {
                someSolutionFound = true;
                base::Cost c = solution->cost(opt_);
                //满足预设就退出，不满足则保留最优结果
                if (opt_->isSatisfied(c))
                {
                    fullyOptimized = true;
                    bestSolution = solution;
                    bestCost_ = c;
                    break;
                }
                else
                {
                    if (opt_->isCostBetterThan(c, bestCost_))
                    {
                        bestSolution = solution;
                        bestCost_ = c;
                        break;  // test cancelled *
                    }
                }
            }
        }
    }

    si_->freeState(workState);

    if (bestSolution)
    {
        base::PlannerSolution psol(bestSolution);
        psol.setPlannerName(getName());
        // if the solution was optimized, we mark it as such
        psol.setOptimized(opt_, bestCost_, fullyOptimized);
        pdef_->addSolutionPath(psol);
    }

    OMPL_INFORM("%s: Created %u states", getName().c_str(), boost::num_vertices(g_) - nrStartStates);

    clock_t end_all = clock();
    if (bestSolution)
    {
        saveLogToFile((end_all - begin_all) * 1000.0 / CLOCKS_PER_SEC, bestCost_.value());
    }
    else
    {
        saveLogToFile(-1, -1);
    }

    return bestSolution ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

void ompl::geometric::LazyPRM::uniteComponents(Vertex a, Vertex b)
{
    unsigned long int componentA = vertexComponentProperty_[a];
    unsigned long int componentB = vertexComponentProperty_[b];
    if (componentA == componentB)
        return;
    if (componentSize_[componentA] > componentSize_[componentB])
    {
        std::swap(componentA, componentB);
        std::swap(a, b);
    }
    markComponent(a, componentB);
}

void ompl::geometric::LazyPRM::markComponent(Vertex v, unsigned long int newComponent)
{
    std::queue<Vertex> q;
    q.push(v);
    while (!q.empty())
    {
        Vertex n = q.front();
        q.pop();
        unsigned long int &component = vertexComponentProperty_[n];
        if (component == newComponent)
            continue;
        if (componentSize_[component] == 1)
            componentSize_.erase(component);
        else
            componentSize_[component]--;
        component = newComponent;
        componentSize_[newComponent]++;
        boost::graph_traits<Graph>::adjacency_iterator nbh, last;
        for (boost::tie(nbh, last) = boost::adjacent_vertices(n, g_); nbh != last; ++nbh)
            q.push(*nbh);
    }
}

//遍历寻找一对已连接的初末姿态
long int ompl::geometric::LazyPRM::solutionComponent(std::pair<std::size_t, std::size_t> *startGoalPair) const
{
    for (std::size_t startIndex = 0; startIndex < startM_.size(); ++startIndex)
    {
        long int startComponent = vertexComponentProperty_[startM_[startIndex]];
        for (std::size_t goalIndex = 0; goalIndex < goalM_.size(); ++goalIndex)
        {
            if (startComponent == (long int)vertexComponentProperty_[goalM_[goalIndex]])
            {
                startGoalPair->first = startIndex;
                startGoalPair->second = goalIndex;
                return startComponent;
            }
        }
    }
    return -1;
}

ompl::base::PathPtr ompl::geometric::LazyPRM::constructSolution(const Vertex &start, const Vertex &goal)
{
    // Need to update the index map here, becuse nodes may have been removed and
    // the numbering will not be 0 .. N-1 otherwise.
    unsigned long int index = 0;
    boost::graph_traits<Graph>::vertex_iterator vi, vend;
    //遍历图中的所有顶点
    for (boost::tie(vi, vend) = boost::vertices(g_); vi != vend; ++vi, ++index)
        indexProperty_[*vi] = index;  //更新索引值

    boost::property_map<Graph, boost::vertex_predecessor_t>::type prev;
    //使用A*算法找路径
    try
    {
        // Consider using a persistent distance_map if it's slow
        boost::astar_search(g_, start, std::bind(&LazyPRM::costHeuristic, this, std::placeholders::_1, goal),
                            boost::predecessor_map(prev)
                                .distance_compare(std::bind(&base::OptimizationObjective::isCostBetterThan, opt_.get(),
                                                            std::placeholders::_1, std::placeholders::_2))
                                .distance_combine(std::bind(&base::OptimizationObjective::combineCosts, opt_.get(),
                                                            std::placeholders::_1, std::placeholders::_2))
                                .distance_inf(opt_->infiniteCost())
                                .distance_zero(opt_->identityCost())
                                .visitor(AStarGoalVisitor<Vertex>(goal)));
    }
    catch (AStarFoundGoal &)
    {
    }
    if (prev[goal] == goal)
        throw Exception(name_, "Could not find solution path");  //?

    // First, get the solution states without copying them, and check them for validity.
    // We do all the node validity checks for the vertices, as this may remove a larger
    // part of the graph (compared to removing an edge).
    std::vector<const base::State *> states(1, stateProperty_[goal]);
    std::set<Vertex> milestonesToRemove;
    //向前递归
    for (Vertex pos = prev[goal]; prev[pos] != pos; pos = prev[pos])
    {
        const base::State *st = stateProperty_[pos];
        unsigned int &vd = vertexValidityProperty_[pos];
        if ((vd & VALIDITY_TRUE) == 0)  // if vd == 0 (false) //检查节点本身是否合法，做标记
            if (si_->isValid(st))
                vd |= VALIDITY_TRUE;
        if ((vd & VALIDITY_TRUE) == 0)
            milestonesToRemove.insert(pos);  //不合法就加入移除集合
        if (milestonesToRemove.empty())      //存储从后到前的最长全部合法路径点集
            states.push_back(st);
    }

    // We remove *all* invalid vertices. This is not entirely as described in the original LazyPRM
    // paper, as the paper suggest removing the first vertex only, and then recomputing the
    // shortest path. However, the paper says the focus is on efficient vertex & edge removal,
    // rather than collision checking, so this modification is in the spirit of the paper.

    //此路径不合法，移除不合法点集，return
    if (!milestonesToRemove.empty())
    {
        unsigned long int comp = vertexComponentProperty_[start];
        // Remember the current neighbors.
        std::set<Vertex> neighbors;
        for (std::set<Vertex>::iterator it = milestonesToRemove.begin(); it != milestonesToRemove.end(); ++it)
        {
            boost::graph_traits<Graph>::adjacency_iterator nbh, last;
            for (boost::tie(nbh, last) = boost::adjacent_vertices(*it, g_); nbh != last; ++nbh)
                if (milestonesToRemove.find(*nbh) == milestonesToRemove.end())
                    neighbors.insert(*nbh);
            // Remove vertex from nearest neighbors data structure.
            nn_->remove(*it);
            // Free vertex state.
            si_->freeState(stateProperty_[*it]);
            // Remove all edges.
            boost::clear_vertex(*it, g_);
            // Remove the vertex.
            boost::remove_vertex(*it, g_);
        }
        // Update the connected component ID for neighbors.
        for (std::set<Vertex>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            if (comp == vertexComponentProperty_[*it])
            {
                unsigned long int newComponent = componentCount_++;
                componentSize_[newComponent] = 0;
                markComponent(*it, newComponent);
            }
        }
        return base::PathPtr();
    }

    // start is checked for validity already
    states.push_back(stateProperty_[start]);

    // Check the edges too, if the vertices were valid. Remove the first invalid edge only.
    //检查边，不合法则去掉第一个不合法的边，并返回
    std::vector<const base::State *>::const_iterator prevState = states.begin(), state = prevState + 1;
    Vertex prevVertex = goal, pos = prev[goal];
    do
    {
        Edge e = boost::lookup_edge(pos, prevVertex, g_).first;
        unsigned int &evd = edgeValidityProperty_[e];
        if ((evd & VALIDITY_TRUE) == 0)
        {
            if (si_->checkMotion(*state, *prevState))
                evd |= VALIDITY_TRUE;
        }
        if ((evd & VALIDITY_TRUE) == 0)
        {
            boost::remove_edge(e, g_);
            unsigned long int newComponent = componentCount_++;
            componentSize_[newComponent] = 0;
            markComponent(pos, newComponent);
            return base::PathPtr();
        }
        prevState = state;
        ++state;
        prevVertex = pos;
        pos = prev[pos];
    } while (prevVertex != pos);
    //规划成功，返回规划结果
    PathGeometric *p = new PathGeometric(si_);
    for (std::vector<const base::State *>::const_reverse_iterator st = states.rbegin(); st != states.rend(); ++st)
        p->append(*st);
    return base::PathPtr(p);
}

ompl::base::Cost ompl::geometric::LazyPRM::costHeuristic(Vertex u, Vertex v) const
{
    return opt_->motionCostHeuristic(stateProperty_[u], stateProperty_[v]);
}

void ompl::geometric::LazyPRM::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    // Explicitly add start and goal states. Tag all states known to be valid as 1.
    // Unchecked states are tagged as 0.
    for (size_t i = 0; i < startM_.size(); ++i)
        data.addStartVertex(base::PlannerDataVertex(stateProperty_[startM_[i]], 1));

    for (size_t i = 0; i < goalM_.size(); ++i)
        data.addGoalVertex(base::PlannerDataVertex(stateProperty_[goalM_[i]], 1));

    // Adding edges and all other vertices simultaneously
    foreach (const Edge e, boost::edges(g_))
    {
        const Vertex v1 = boost::source(e, g_);
        const Vertex v2 = boost::target(e, g_);
        data.addEdge(base::PlannerDataVertex(stateProperty_[v1]), base::PlannerDataVertex(stateProperty_[v2]));

        // Add the reverse edge, since we're constructing an undirected roadmap
        data.addEdge(base::PlannerDataVertex(stateProperty_[v2]), base::PlannerDataVertex(stateProperty_[v1]));

        // Add tags for the newly added vertices
        data.tagState(stateProperty_[v1], (vertexValidityProperty_[v1] & VALIDITY_TRUE) == 0 ? 0 : 1);
        data.tagState(stateProperty_[v2], (vertexValidityProperty_[v2] & VALIDITY_TRUE) == 0 ? 0 : 1);
    }
}

bool ompl::geometric::LazyPRM::saveLogToFile(double time, double cost)
{
    std::string home_path = getenv("HOME");
    std::string file_name_path = "/tmp/rm_name";
    std::string filename;
    std::fstream namefin(file_name_path, std::ios::in);
    if (!namefin.is_open())
    {
        OMPL_ERROR("unable to open file %s", (file_name_path).c_str());
    }
    namefin >> filename;
    namefin.close();
    printf("writing to file...");
    std::string save_path_full = "/mgn_data/test_log/LazyPRMlog.txt";
    std::fstream fout(home_path + save_path_full, std::ios::app);
    if (!fout.is_open())
    {
        std::cerr << "unable to open file " << home_path + save_path_full << std::endl;
    }
    fout << filename << "  time " << time << "  cost " << cost;
    std::cout << filename << "  time " << time << "  cost " << cost;
    fout << std::endl;
    fout.close();
    std::cout << std::endl;
    return true;
}

int ompl::geometric::LazyPRM::addGeneratdMilestones()
{
    std::string home_path = getenv("HOME");
    double pre_map_time_max = 50;
    clock_t begin_pre_map = clock();
    clock_t end_pre_map = clock();

    OMPL_DEBUG("addGeneratdMilestones");
    base::State *state = si_->allocState();
    std::string file_name_path = "/tmp/rm_name";
    std::string filename;
    std::fstream namefin(file_name_path, std::ios::in);
    if (!namefin.is_open())
    {
        OMPL_ERROR("unable to open file %s", file_name_path.c_str());
        return false;
    }
    namefin >> filename;
    namefin.close();

    OMPL_INFORM("filename %s", filename.c_str());
    std::string filenamefullpath = "/mgn_data/randmat6d.txt";
    std::fstream fin(home_path + filenamefullpath, std::ios::in);
    //    fout.open(filename_.data(),ios::in|ios::out);
    //    fout.open("filename_toFile.txt",ios::in|ios::out);
    if (!fin.is_open())
    {
        OMPL_ERROR("unable to open file %s check path", filenamefullpath.c_str());
        return false;
    }
    char buffer[256];
    int vertexNum = 0;

    int flag_max = 100;  //按一定比例加入生成点和随机点
    int flag = flag_max;
    int cnt_max = 512;
    while (!fin.eof() && (cnt_max-- > 0) &&
           ((end_pre_map - begin_pre_map) * 1000.0 / CLOCKS_PER_SEC) < pre_map_time_max)
    {
        flag--;
        if (flag < 1)
        {
            if (!sampler_)
            {
                sampler_ = si_->allocStateSampler();
            }
            sampler_->sampleUniform(state);
            flag = flag_max;
            //            OMPL_INFORM("sample");
        }
        else
        {
            double *val = static_cast<ompl::base::RealVectorStateSpace::StateType *>(state)->values;
            //        double point[3];
            fin.getline(buffer, 100);
            sscanf(buffer, "%lf %lf %lf %lf %lf %lf\n", &val[0], &val[1], &val[2], &val[3], &val[4], &val[5]);
            //            printf("%lf %lf %lf\n", val[0], val[1], val[2]);
            if (fabs(val[0]) > 5.0 || fabs(val[1]) > 5.0 || fabs(val[2]) > 5.0)
            {
                continue;
            }
            // OMPL_INFORM("load");
        }
        //        for(int i=0;i<3;i++){
        //
        //        }
        //        OMPL_INFORM("%f %f %f", val[0], val[1], val[2]);
        // OMPL_INFORM("sample");
        /*Vertex addedVertex =*/
        if (si_->isValid(state))
        {
            addMilestone(si_->cloneState(state));
            // OMPL_INFORM("add");
            vertexNum++;
        }
        end_pre_map = clock();
        //        fout << val[0] << "  " << val[1] << "  " << val[2] << endl;
    }
    si_->freeState(state);
    OMPL_INFORM("Vertex number loaded: %d", vertexNum);
    OMPL_INFORM("Vertex pre-map time: %lf", (end_pre_map - begin_pre_map) * 1000.0 / CLOCKS_PER_SEC);
    fin.close();

    return true;
}