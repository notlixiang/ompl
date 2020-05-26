/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011, Rice University
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
 *   * Neither the name of Rice University nor the names of its
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

/* Author: Ioan Sucan, James D. Marble, Ryan Luna */

#include "ompl/geometric/planners/prm/ValidStateGen.h"
#include "ompl/geometric/planners/prm/ConnectionStrategy.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/datastructures/PDF.h"
#include "ompl/tools/config/SelfConfig.h"
#include "ompl/tools/config/MagicConstants.h"
#include <boost/graph/astar_search.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/foreach.hpp>
#include <thread>

#include "GoalVisitor.hpp"

#define foreach BOOST_FOREACH
#define MAX_TRY 1000

namespace ompl
{
    namespace magic
    {
        /** \brief The number of steps to take for a random bounce
                motion generated as part of the expansion step of ValidStateGen. */
        static const unsigned int MAX_RANDOM_BOUNCE_STEPS = 5;

        /** \brief The time in seconds for a single roadmap building operation (dt)*/
        static const double ROADMAP_BUILD_TIME = 0.2;

        /** \brief The number of nearest neighbors to consider by
                default in the construction of the ValidStateGen roadmap */
        static const unsigned int DEFAULT_NEAREST_NEIGHBORS = 10;
    }  // namespace magic
}  // namespace ompl

// double ompl::base::RealVectorStateSpace::distance(const State *state1, const State *state2) const
//{
//    double dist = 0.0;
//    const double *s1 = static_cast<const StateType*>(state1)->values;
//    const double *s2 = static_cast<const StateType*>(state2)->values;
//
//    for (unsigned int i = 0 ; i < dimension_ ; ++i)
//    {
//        double diff = (*s1++) - (*s2++);
//        dist += diff * diff;
//    }
//    return sqrt(dist);
//}

ompl::geometric::ValidStateGen::ValidStateGen(const base::SpaceInformationPtr &si, bool starStrategy)
  : base::Planner(si, "ValidStateGen")
  , starStrategy_(starStrategy)
  , stateProperty_(boost::get(vertex_state_t(), g_))
  , totalConnectionAttemptsProperty_(boost::get(vertex_total_connection_attempts_t(), g_))
  , successfulConnectionAttemptsProperty_(boost::get(vertex_successful_connection_attempts_t(), g_))
  , weightProperty_(boost::get(boost::edge_weight, g_))
  , disjointSets_(boost::get(boost::vertex_rank, g_), boost::get(boost::vertex_predecessor, g_))
  , userSetConnectionStrategy_(false)
  , addedNewSolution_(false)
  , iterations_(0)
  , bestCost_(std::numeric_limits<double>::quiet_NaN())
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.approximateSolutions = false;
    specs_.optimizingPaths = true;
    specs_.multithreaded = true;

    if (!starStrategy_)
        Planner::declareParam<unsigned int>("max_nearest_neighbors", this, &ValidStateGen::setMaxNearestNeighbors,
                                            std::string("8:1000"));

    addPlannerProgressProperty("iterations INTEGER", std::bind(&ValidStateGen::getIterationCount, this));
    addPlannerProgressProperty("best cost REAL", std::bind(&ValidStateGen::getBestCost, this));
    addPlannerProgressProperty("milestone count INTEGER", std::bind(&ValidStateGen::getMilestoneCountString, this));
    addPlannerProgressProperty("edge count INTEGER", std::bind(&ValidStateGen::getEdgeCountString, this));
}

ompl::geometric::ValidStateGen::~ValidStateGen()
{
    freeMemory();
}

void ompl::geometric::ValidStateGen::setup()
{
    Planner::setup();
    if (!nn_)
    {
        specs_.multithreaded = false;  // temporarily set to false since nn_ is used only in single thread
        nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
        specs_.multithreaded = true;
        nn_->setDistanceFunction(
            std::bind(&ValidStateGen::distanceFunction, this, std::placeholders::_1, std::placeholders::_2));
    }
    if (!connectionStrategy_)
    {
        if (starStrategy_)
            connectionStrategy_ =
                KStarStrategy<Vertex>(std::bind(&ValidStateGen::milestoneCount, this), nn_, si_->getStateDimension());
        else
            connectionStrategy_ = KStrategy<Vertex>(magic::DEFAULT_NEAREST_NEIGHBORS, nn_);
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
        {
            OMPL_INFORM("hasOptimizationObjective");
            opt_ = pdef_->getOptimizationObjective();
        }

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
}

void ompl::geometric::ValidStateGen::setMaxNearestNeighbors(unsigned int k)
{
    if (starStrategy_)
        throw Exception("Cannot set the maximum nearest neighbors for " + getName());
    if (!nn_)
    {
        specs_.multithreaded = false;  // temporarily set to false since nn_ is used only in single thread
        nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
        specs_.multithreaded = true;
        nn_->setDistanceFunction(
            std::bind(&ValidStateGen::distanceFunction, this, std::placeholders::_1, std::placeholders::_2));
    }
    if (!userSetConnectionStrategy_)
        connectionStrategy_ = ConnectionStrategy();
    if (isSetup())
        setup();
}

void ompl::geometric::ValidStateGen::setProblemDefinition(const base::ProblemDefinitionPtr &pdef)
{
    Planner::setProblemDefinition(pdef);
    clearQuery();
}

void ompl::geometric::ValidStateGen::clearQuery()
{
    startM_.clear();
    goalM_.clear();
    pis_.restart();
}

void ompl::geometric::ValidStateGen::clear()
{
    Planner::clear();
    sampler_.reset();
    simpleSampler_.reset();
    freeMemory();
    if (nn_)
        nn_->clear();
    clearQuery();

    iterations_ = 0;
    bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
}

void ompl::geometric::ValidStateGen::freeMemory()
{
    foreach (Vertex v, boost::vertices(g_))
        si_->freeState(stateProperty_[v]);
    g_.clear();
}

void ompl::geometric::ValidStateGen::expandRoadmap(double expandTime)
{
    expandRoadmap(base::timedPlannerTerminationCondition(expandTime));
}

void ompl::geometric::ValidStateGen::expandRoadmap(const base::PlannerTerminationCondition &ptc)
{
    if (!simpleSampler_)
        simpleSampler_ = si_->allocStateSampler();

    std::vector<base::State *> states(magic::MAX_RANDOM_BOUNCE_STEPS);
    si_->allocStates(states);
    expandRoadmap(ptc, states);
    si_->freeStates(states);
}

void ompl::geometric::ValidStateGen::expandRoadmap(const base::PlannerTerminationCondition &ptc,
                                                   std::vector<base::State *> &workStates)
{
    // construct a probability distribution over the vertices in the roadmap
    // as indicated in
    //  "Probabilistic Roadmaps for Path Planning in High-Dimensional Configuration Spaces"
    //        Lydia E. Kavraki, Petr Svestka, Jean-Claude Latombe, and Mark H. Overmars

    PDF<Vertex> pdf;
    foreach (Vertex v, boost::vertices(g_))
    {
        const unsigned long int t = totalConnectionAttemptsProperty_[v];
        pdf.add(v, (double)(t - successfulConnectionAttemptsProperty_[v]) / (double)t);
        //失败连接次数占总次数的比例，越高说明此处环境越复杂
    }

    if (pdf.empty())
        return;  //在足够的growRoadmap才能使用expandRoadmap

    while (ptc == false)
    {
        iterations_++;
        Vertex v = pdf.sample(rng_.uniform01());  //使用PDF采样
        unsigned int s =
            si_->randomBounceMotion(simpleSampler_, stateProperty_[v], workStates.size(), workStates, false);
        //随机尝试连接，返回试出的合法路径
        // workStates用来存储结果,已分配空间不需要重新分配
        if (s > 0)
        {
            s--;
            OMPL_INFORM("expandRoadmap");
            Vertex last = addMilestone(si_->cloneState(workStates[s]));

            graphMutex_.lock();
            for (unsigned int i = 0; i < s; ++i)
            {
                // add the vertex along the bouncing motion
                Vertex m = boost::add_vertex(g_);
                stateProperty_[m] = si_->cloneState(workStates[i]);
                totalConnectionAttemptsProperty_[m] = 1;
                successfulConnectionAttemptsProperty_[m] = 0;
                disjointSets_.make_set(m);

                // add the edge to the parent vertex
                const base::Cost weight = opt_->motionCost(stateProperty_[v], stateProperty_[m]);
                const Graph::edge_property_type properties(weight);
                boost::add_edge(v, m, properties, g_);
                uniteComponents(v, m);

                // add the vertex to the nearest neighbors data structure
                nn_->add(m);
                v = m;
            }

            // if there are intermediary states or the milestone has not been connected to the initially sampled vertex,
            // we add an edge
            if (s > 0 || !sameComponent(v, last))
            {
                // add the edge to the parent vertex
                const base::Cost weight = opt_->motionCost(stateProperty_[v], stateProperty_[last]);
                const Graph::edge_property_type properties(weight);
                boost::add_edge(v, last, properties, g_);
                uniteComponents(v, last);
            }
            graphMutex_.unlock();
        }
    }
}

void ompl::geometric::ValidStateGen::growRoadmap(double growTime)
{
    OMPL_INFORM("growRoadmap(double growTime)");
    growRoadmap(base::timedPlannerTerminationCondition(growTime));
}

void ompl::geometric::ValidStateGen::growRoadmap(const base::PlannerTerminationCondition &ptc)
{
    OMPL_INFORM("growRoadmap(const base::PlannerTerminationCondition &ptc)");
    if (!isSetup())
        setup();
    if (!sampler_)
        sampler_ = si_->allocValidStateSampler();

    base::State *workState = si_->allocState();

    growRoadmap(ptc, workState);
    si_->freeState(workState);
}

void ompl::geometric::ValidStateGen::growRoadmap(const base::PlannerTerminationCondition &ptc, base::State *workState)
{
    /* grow roadmap in the regular fashion -- sample valid states, add them to the roadmap, add valid connections */
    while (ptc == false)
    {
        iterations_++;
        // search for a valid state
        bool found = false;
        while (!found && ptc == false)
        {
            unsigned int attempts = 0;
            do
            {
                found = sampler_->sample(workState);
                attempts++;
            } while (attempts < magic::FIND_VALID_STATE_ATTEMPTS_WITHOUT_TERMINATION_CHECK && !found);
        }
        // add it as a milestone
        if (found)
        {
            //            OMPL_INFORM("growRoadmap");
            addMilestone(si_->cloneState(workState));
        }
    }
}

void ompl::geometric::ValidStateGen::checkForSolution(const base::PlannerTerminationCondition &ptc,
                                                      base::PathPtr &solution)
{
    base::GoalSampleableRegion *goal = static_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());
    while (!ptc && !addedNewSolution_)
    {
        // Check for any new goal states
        if (goal->maxSampleCount() > goalM_.size())
        {
            const base::State *st = pis_.nextGoal();
            if (st)
            {
                OMPL_INFORM("checkForSolution");
                goalM_.push_back(addMilestone(si_->cloneState(st)));
            }
        }

        // Check for a solution
        addedNewSolution_ = maybeConstructSolution(startM_, goalM_, solution);
        // Sleep for 1ms
        if (!addedNewSolution_)
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

bool ompl::geometric::ValidStateGen::maybeConstructSolution(const std::vector<Vertex> &starts,
                                                            const std::vector<Vertex> &goals, base::PathPtr &solution)
{
    base::Goal *g = pdef_->getGoal().get();
    base::Cost sol_cost(opt_->infiniteCost());
    foreach (Vertex start, starts)
    {
        foreach (Vertex goal, goals)
        {
            // we lock because the connected components algorithm is incremental and may change disjointSets_
            graphMutex_.lock();
            bool same_component = sameComponent(start, goal);  //检查是否在同一连通子图上
            graphMutex_.unlock();

            if (same_component && g->isStartGoalPairValid(stateProperty_[goal], stateProperty_[start]))
            {
                base::PathPtr p = constructSolution(start, goal);
                if (p)
                {
                    base::Cost pathCost = p->cost(opt_);
                    if (opt_->isCostBetterThan(pathCost, bestCost_))
                        bestCost_ = pathCost;
                    // Check if optimization objective is satisfied
                    if (opt_->isSatisfied(pathCost))
                    {
                        solution = p;
                        return true;
                    }
                    else if (opt_->isCostBetterThan(pathCost, sol_cost))
                    {
                        solution = p;
                        sol_cost = pathCost;
                    }
                }
            }
        }
    }

    return false;
}

bool ompl::geometric::ValidStateGen::addedNewSolution() const
{
    return addedNewSolution_;
}

ompl::base::PlannerStatus ompl::geometric::ValidStateGen::solve(const base::PlannerTerminationCondition &ptc)
{
    OMPL_INFORM("Start generation ...");
    checkValidity();
    base::GoalSampleableRegion *goal = dynamic_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());

    if (!goal)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }

    bool sampable = testSampable(100);
    if (!sampable)
    {
        getVertexData();
        saveTryTimes();
        return base::PlannerStatus::TIMEOUT;
    }

    // Add the valid start states as milestones
    //    while (const base::State *st = pis_.nextStart())
    //        startM_.push_back(addMilestone(si_->cloneState(st)));
    //
    //    if (startM_.size() == 0) {
    //        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
    //        return base::PlannerStatus::INVALID_START;
    //    }

    ////不需要检查初末姿态合法性

    //    if (!goal->couldSample()) {
    //        OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
    //        return base::PlannerStatus::INVALID_GOAL;
    //    }
    //
    //    // Ensure there is at least one valid goal state
    //    if (goal->maxSampleCount() > goalM_.size() || goalM_.empty()) {
    //        const base::State *st = goalM_.empty() ? pis_.nextGoal(ptc) : pis_.nextGoal();
    //        if (st)
    //            goalM_.push_back(addMilestone(si_->cloneState(st)));
    //
    //        if (goalM_.empty()) {
    //            OMPL_ERROR("%s: Unable to find any valid goal states", getName().c_str());
    //            return base::PlannerStatus::INVALID_GOAL;
    //        }
    //    }
    //
    unsigned long int nrStartStates = boost::num_vertices(g_);
    OMPL_INFORM("%s: Starting planning with %lu states already in datastructure", getName().c_str(), nrStartStates);
    //以上都在初始化

    // Reset addedNewSolution_ member and create solution checking thread
    //    addedNewSolution_ = false;
    //    base::PathPtr sol;

    // construct new planner termination condition that fires when the given ptc is true, or a solution is found
    //    base::PlannerTerminationCondition ptcOrSolutionFound =
    //            base::plannerOrTerminationCondition(ptc, base::PlannerTerminationCondition(
    //                    std::bind(&ValidStateGen::addedNewSolution, this)));

    try_times.clear();

    std::thread slnThread1(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread2(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread3(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread4(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread5(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread6(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread7(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));
    //    std::thread slnThread8(std::bind(&ValidStateGen::addValidVertexThread, this, ptc, MAX_TRY));

    slnThread1.join();
    //    slnThread2.join();
    //    slnThread3.join();
    //    slnThread4.join();
    //    slnThread5.join();
    //    slnThread6.join();
    //    slnThread7.join();
    //    slnThread8.join();

    //    int iter_try_add_main = 0;
    //    while (ptc() == false) {
    //        int try_num = addValidVertex(ptc, MAX_TRY);
    //        OMPL_INFORM("try_num : %5d", try_num);
    //        OMPL_INFORM("iter_try_add_main : %5d", iter_try_add_main++);
    //        try_times.push_back(try_num);
    //    }
    //    while (ptc() == false) { ;
    //    }

    getVertexData();
    saveTryTimes();

    //    // Ensure slnThread is ceased before exiting solve

    //    OMPL_INFORM("%s: Created %u states", getName().c_str(), boost::num_vertices(g_) - nrStartStates);
    //
    //    if (sol) {
    //        base::PlannerSolution psol(sol);
    //        psol.setPlannerName(getName());
    //        // if the solution was optimized, we mark it as such
    //        psol.setOptimized(opt_, bestCost_, addedNewSolution());
    //        pdef_->addSolutionPath(psol);
    //    }

    return base::PlannerStatus::TIMEOUT;
}

void ompl::geometric::ValidStateGen::constructRoadmap(const base::PlannerTerminationCondition &ptc)
{
    if (!isSetup())
        setup();
    if (!sampler_)
        sampler_ = si_->allocValidStateSampler();
    if (!simpleSampler_)
        simpleSampler_ = si_->allocStateSampler();

    std::vector<base::State *> xstates(magic::MAX_RANDOM_BOUNCE_STEPS);
    si_->allocStates(xstates);
    bool grow = true;

    bestCost_ = opt_->infiniteCost();
    while (ptc() == false)
    {
        // maintain a 2:1 ratio for growing/expansion of roadmap
        // call growRoadmap() twice as long for every call of expandRoadmap()
        if (grow)
        {
            OMPL_INFORM("constructRoadmap");
            growRoadmap(base::plannerOrTerminationCondition(
                            ptc, base::timedPlannerTerminationCondition(2.0 * magic::ROADMAP_BUILD_TIME)),
                        xstates[0]);
        }
        else
            expandRoadmap(base::plannerOrTerminationCondition(
                              ptc, base::timedPlannerTerminationCondition(magic::ROADMAP_BUILD_TIME)),
                          xstates);
        grow = !grow;
    }

    si_->freeStates(xstates);
}

ompl::geometric::ValidStateGen::Vertex ompl::geometric::ValidStateGen::addMilestone(base::State *state)
{
    std::lock_guard<std::mutex> _(graphMutex_);

    double *val = static_cast<ompl::base::RealVectorStateSpace::StateType *>(state)->values;
    //    OMPL_INFORM("%f %f %f", val[0], val[1], val[2]);
    OMPL_INFORM("addMilestone");

    Vertex m = boost::add_vertex(g_);
    stateProperty_[m] = state;
    totalConnectionAttemptsProperty_[m] = 1;
    successfulConnectionAttemptsProperty_[m] = 0;

    // Initialize to its own (dis)connected component.
    disjointSets_.make_set(m);

    // Which milestones will we attempt to connect to?
    const std::vector<Vertex> &neighbors = connectionStrategy_(m);

    foreach (Vertex n, neighbors)  // for each neighbour
        if (connectionFilter_(n, m))
        {  // always true?
            totalConnectionAttemptsProperty_[m]++;
            totalConnectionAttemptsProperty_[n]++;
            if (si_->checkMotion(stateProperty_[n], stateProperty_[m]))
            {
                successfulConnectionAttemptsProperty_[m]++;
                successfulConnectionAttemptsProperty_[n]++;
                const base::Cost weight = opt_->motionCost(stateProperty_[n], stateProperty_[m]);
                const Graph::edge_property_type properties(weight);
                boost::add_edge(n, m, properties, g_);
                uniteComponents(n, m);
            }
        }

    nn_->add(m);

    return m;
}

void ompl::geometric::ValidStateGen::uniteComponents(Vertex m1, Vertex m2)
{
    disjointSets_.union_set(m1, m2);
}

bool ompl::geometric::ValidStateGen::sameComponent(Vertex m1, Vertex m2)
{
    return boost::same_component(m1, m2, disjointSets_);
}

ompl::base::PathPtr ompl::geometric::ValidStateGen::constructSolution(const Vertex &start, const Vertex &goal)
{
    std::lock_guard<std::mutex> _(graphMutex_);
    boost::vector_property_map<Vertex> prev(boost::num_vertices(g_));

    try
    {
        // Consider using a persistent distance_map if it's slow
        boost::astar_search(g_, start, std::bind(&ValidStateGen::costHeuristic, this, std::placeholders::_1, goal),
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
        throw Exception(name_, "Could not find solution path");

    PathGeometric *p = new PathGeometric(si_);
    for (Vertex pos = goal; prev[pos] != pos; pos = prev[pos])
        p->append(stateProperty_[pos]);
    p->append(stateProperty_[start]);
    p->reverse();

    return base::PathPtr(p);
}

void ompl::geometric::ValidStateGen::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    // Explicitly add start and goal states:
    for (size_t i = 0; i < startM_.size(); ++i)
        data.addStartVertex(base::PlannerDataVertex(
            stateProperty_[startM_[i]], const_cast<ValidStateGen *>(this)->disjointSets_.find_set(startM_[i])));

    for (size_t i = 0; i < goalM_.size(); ++i)
        data.addGoalVertex(base::PlannerDataVertex(
            stateProperty_[goalM_[i]], const_cast<ValidStateGen *>(this)->disjointSets_.find_set(goalM_[i])));

    // Adding edges and all other vertices simultaneously
    foreach (const Edge e, boost::edges(g_))
    {
        const Vertex v1 = boost::source(e, g_);
        const Vertex v2 = boost::target(e, g_);
        data.addEdge(base::PlannerDataVertex(stateProperty_[v1]), base::PlannerDataVertex(stateProperty_[v2]));

        // Add the reverse edge, since we're constructing an undirected roadmap
        data.addEdge(base::PlannerDataVertex(stateProperty_[v2]), base::PlannerDataVertex(stateProperty_[v1]));

        // Add tags for the newly added vertices
        data.tagState(stateProperty_[v1], const_cast<ValidStateGen *>(this)->disjointSets_.find_set(v1));
        data.tagState(stateProperty_[v2], const_cast<ValidStateGen *>(this)->disjointSets_.find_set(v2));
    }
}

ompl::base::Cost ompl::geometric::ValidStateGen::costHeuristic(Vertex u, Vertex v) const
{
    return opt_->motionCostHeuristic(stateProperty_[u], stateProperty_[v]);
}

bool ompl::geometric::ValidStateGen::tryAddValidAndSparse(Vertex m, base::State *state)
{
    clock_t begin = clock();
    clock_t end = clock();
    OMPL_INFORM("tryAddValidAndSparse");
    std::lock_guard<std::mutex> _(graphMutex_);

    double *val = static_cast<ompl::base::RealVectorStateSpace::StateType *>(state)->values;
    OMPL_INFORM("%f %f %f %f %f %f", val[0], val[1], val[2], val[3], val[4], val[5]);

    //    Vertex m = boost::add_vertex(g_);
    stateProperty_[m] = state;
    totalConnectionAttemptsProperty_[m] = 1;
    successfulConnectionAttemptsProperty_[m] = 0;

    // Initialize to its own (dis)connected component.

    //    disjointSets_.make_set(m);//add this to upper function

    // Which milestones will we attempt to connect to?
    const std::vector<Vertex> &neighbors = connectionStrategy_(m);
    std::vector<Vertex> neighborsConnectable;
    neighborsConnectable.clear();
    bool flag = false;
    end = clock();
    OMPL_INFORM("setup time: %f ms", (end - begin) * 1000.0 / CLOCKS_PER_SEC);
    begin = clock();
    BOOST_FOREACH (Vertex n, neighbors)  // for each neighbour
    {
        if (connectionFilter_(n, m))
        {  // always true?
            //                        totalConnectionAttemptsProperty_[m]++;
            //                        totalConnectionAttemptsProperty_[n]++;
            if (si_->checkMotion(stateProperty_[n], stateProperty_[m]))
            {
                neighborsConnectable.push_back(n);
                //                            successfulConnectionAttemptsProperty_[m]++;
                //                            successfulConnectionAttemptsProperty_[n]++;
                //                            const base::Cost weight = opt_->motionCost(stateProperty_[n],
                //                            stateProperty_[m]); const Graph::edge_property_type properties(weight);
                //                            boost::add_edge(n, m, properties, g_);
                //                            uniteComponents(n, m);
            }
        }
    }

    end = clock();
    OMPL_INFORM("neighbor connectable check time: %f ms", (end - begin) * 1000.0 / CLOCKS_PER_SEC);
    if (neighborsConnectable.size() < 2)
    {  //可连数过小，场景复杂，加点
        flag = true;
        //        OMPL_INFORM("neighborsConnectable.size() %ld",neighborsConnectable.size());
    }
    else
    {  //若可连的点不全在同一子图上，加点
        begin = clock();
        bool flag_existDifferentComponent = false;
        BOOST_FOREACH (Vertex n1, neighborsConnectable)
        {
            BOOST_FOREACH (Vertex n2, neighborsConnectable)
            {
                //                                        graphMutex_.lock();
                bool same_component = sameComponent(n1, n2);  //检查是否在同一连通子图上
                //                                        graphMutex_.unlock();
                if (same_component == false)
                {
                    flag_existDifferentComponent = true;
                }
                if (flag_existDifferentComponent == true)
                {
                    break;
                }
            }
            if (flag_existDifferentComponent == true)
            {
                break;
            }
        }
        end = clock();
        OMPL_INFORM("neighbor sameComponent check time: %f ms", (end - begin) * 1000.0 / CLOCKS_PER_SEC);
        if (flag_existDifferentComponent == true)
        {
            flag = true;
            //            OMPL_INFORM("neighborsConnectable.size() %ld",neighborsConnectable.size());
        }
        else
        {  // n点相互连接
            begin = clock();
            bool flag_existUnconnectableNeib = false;
            size_t len = neighborsConnectable.size();
            for (size_t n1 = 0; n1 < len; n1++)
            {
                for (size_t n2 = 0; n2 < n1; n2++)
                {
                    bool connectable = si_->checkMotion(stateProperty_[neighborsConnectable[n1]],
                                                        stateProperty_[neighborsConnectable[n2]]);  //检查是否能够直连
                    if (connectable == false)
                    {
                        flag_existUnconnectableNeib = true;
                    }
                    if (flag_existUnconnectableNeib == true)
                    {
                        break;
                    }
                }
                if (flag_existUnconnectableNeib == true)
                {
                    break;
                }
            }
            end = clock();
            OMPL_INFORM("neighbor checkMotion time: %f ms", (end - begin) * 1000.0 / CLOCKS_PER_SEC);
            if (flag_existUnconnectableNeib == true)
            {
                flag = true;
            }
        }
    }

    if (flag == false)
    {
        return false;
        OMPL_INFORM("tryAddValidAndSparse stop,false");
    }
    else
    {
        begin = clock();
        BOOST_FOREACH (Vertex n, neighborsConnectable)  // for each neighbour
        {
            if (connectionFilter_(n, m))
            {  // always true?
                //                        totalConnectionAttemptsProperty_[m]++;
                //                        totalConnectionAttemptsProperty_[n]++;

                successfulConnectionAttemptsProperty_[m]++;
                successfulConnectionAttemptsProperty_[n]++;
                const base::Cost weight = opt_->motionCost(stateProperty_[n], stateProperty_[m]);
                // OMPL_INFORM("weight %f",weight.value());
                //                            std::cout<<"weight: "<<weight<<std::endl;
                const Graph::edge_property_type properties(weight);
                boost::add_edge(n, m, properties, g_);
                uniteComponents(n, m);
            }
        }
        nn_->add(m);
        end = clock();
        OMPL_INFORM("add vertex time: %f ms", (end - begin) * 1000.0 / CLOCKS_PER_SEC);
        OMPL_INFORM("tryAddValidAndSparse stop,true");
        return true;
    }
}

int ompl::geometric::ValidStateGen::addValidVertex(const base::PlannerTerminationCondition &ptc, int max_try)
{
    OMPL_INFORM("addValidVertex");
    if (!isSetup())
        setup();
    if (!sampler_)
        sampler_ = si_->allocValidStateSampler();
    if (!simpleSampler_)
        simpleSampler_ = si_->allocStateSampler();

    base::State *workState = si_->allocState();

    Vertex m = boost::add_vertex(g_);
    int iter = 0;

    //    growRoadmap(ptc, workState);

    while (ptc == false && iter < max_try)
    {
        OMPL_INFORM("addValidVertex iter %d", iter);
        iter++;
        // search for a valid state
        bool found = false;
        while (!found && ptc == false)
        {
            unsigned int attempts = 0;
            do
            {
                //                OMPL_INFORM("attempts %d", attempts);
                found = sampler_->sample(workState);
                attempts++;
            } while (attempts < magic::FIND_VALID_STATE_ATTEMPTS_WITHOUT_TERMINATION_CHECK && !found);
        }
        // add it as a milestone
        if (found)
        {
            OMPL_INFORM("found");
            if (tryAddValidAndSparse(m, si_->cloneState(workState)))
            {
                break;
            }
        }
    }
    si_->freeState(workState);
    OMPL_INFORM("addValidVertex stop");
    return iter;
}

bool ompl::geometric::ValidStateGen::getVertexData()
{
    OMPL_INFORM("getVertexData");
    std::string home_path = getenv("HOME");
    const char *file_name_path = "/mgn_data/name";
    std::string filename;
    std::fstream fin(home_path + file_name_path, std::ios::in);
    if (!fin.is_open())
    {
        OMPL_ERROR("无法打开文件 %s", file_name_path);
        return false;
    }
    fin >> filename;
    fin.close();

    OMPL_INFORM("filename %s", filename.data());
    std::string filenamefullpath = "/home/null/mgn_data/prm/" + filename + ".prm";
    std::fstream fout(filenamefullpath.data(), std::ios::out);
    //    fout.open(filename_.data(),ios::in|ios::out);
    //    fout.open("filename_toFile.txt",ios::in|ios::out);
    if (!fout.is_open())
    {
        OMPL_ERROR("无法打开文件 %s, 请检查文件路径", filenamefullpath);
        return false;
    }
    Graph::vertex_iterator vertexIt, vertexEnd;  // 顶点
    boost::tie(vertexIt, vertexEnd) = boost::vertices(g_);
    for (; vertexIt != vertexEnd; ++vertexIt)
    {
        base::State *state = stateProperty_[*vertexIt];
        double *val = static_cast<ompl::base::RealVectorStateSpace::StateType *>(state)->values;
        // OMPL_INFORM("%f %f %f %f %f %f", val[0], val[1], val[2], val[3], val[4], val[5]);
        fout << val[0] << "  " << val[1] << "  " << val[2] << "  " << val[3] << "  " << val[4] << "  " << val[5]
             << endl;
    }
    OMPL_INFORM("vertex number: %ld", *vertexEnd);
    fout.close();

    return true;
}

bool ompl::geometric::ValidStateGen::saveTryTimes()
{
    OMPL_INFORM("saveTryTimes");

    std::string home_path = getenv("HOME");
    const char *file_name_path = "/mgn_data/name";
    std::string filename;
    std::fstream fin(home_path + file_name_path, std::ios::in);
    if (!fin.is_open())
    {
        OMPL_ERROR("无法打开文件 %s", file_name_path);
        return false;
    }
    fin >> filename;
    fin.close();

    OMPL_INFORM("filename %s", filename.data());
    std::string filenamefullpath = "/home/null/mgn_data/TryTimes/" + filename + ".pts";
    std::fstream fout(filenamefullpath.data(), std::ios::out);
    //    fout.open(filename_.data(),ios::in|ios::out);
    //    fout.open("filename_toFile.txt",ios::in|ios::out);
    if (!fout.is_open())
    {
        OMPL_ERROR("无法打开文件 %s, 请检查文件路径", filenamefullpath);
        return false;
    }
    size_t len = try_times.size();
    if (len > 0)
    {
        for (int i = 0; i < len; i++)
        {
            fout << try_times[i] << endl;
        }
    }
    fout.close();

    return true;
}

// void ompl::geometric::ValidStateGen::saveCloudToFile(std::string path){}

void ompl::geometric::ValidStateGen::addValidVertexThread(const base::PlannerTerminationCondition &ptc, int max_try)
{
    int iter_try_add_main = 0;
    while (ptc() == false)
    {
        int try_num = addValidVertex(ptc, MAX_TRY);
        OMPL_INFORM("try_num : %5d", try_num);
        OMPL_INFORM("iter_try_add_main : %5d", iter_try_add_main++);
        if (try_num == -1)
        {
            break;
        }
        try_times.push_back(try_num);
    }
}

bool ompl::geometric::ValidStateGen::testSampable(int max_try)
{
    if (!isSetup())
        setup();
    if (!sampler_)
        sampler_ = si_->allocValidStateSampler();
    if (!simpleSampler_)
        simpleSampler_ = si_->allocStateSampler();

    base::State *workState = si_->allocState();

    int iter = 0;

    //    growRoadmap(ptc, workState);

    bool found = false;
    while (iter < max_try)
    {
        iter++;
        // search for a valid state

        found = sampler_->sample(workState);

        // add it as a milestone
        if (found)
        {
            OMPL_INFORM("sampable");
            si_->freeState(workState);
            return true;
        }
    }
    OMPL_INFORM("not sampable");
    si_->freeState(workState);
    return false;
}