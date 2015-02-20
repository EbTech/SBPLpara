/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Pennsylvania nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __EPASEPlanner_H_
#define __EPASEPlanner_H_

#include <cstdio>
#include <ctime>
#include <vector>
#include <map>
#include <sbpl/planners/planner.h>
#include <sbpl/utils/mdp.h>

#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>

//---configuration----
//control of EPS
//initial suboptimality bound (cost solution <= cost(eps*cost optimal solution)
#define EPASE_DEFAULT_INITIAL_EPS	    5.0
//as planning time exist, EPASE* decreases epsilon bound
#define EPASE_DECREASE_EPS    0.2
//final epsilon bound
#define EPASE_FINAL_EPS	    1.0
//---------------------
#define EPASE_INCONS_LIST_ID 0

#define EPASE_NUM_THREADS 4

class CList;
class DiscreteSpaceInformation;
class MDPConfig;
class StateChangeQuery;

//-------------------------------------------------------------

/**
 * \brief state structure used in EPASE* search tree
 */
typedef class EPASESearchStateData : public AbstractSearchState
{
public:
    /**
     * \brief the MDP state itself
     */
    CMDPSTATE* MDPstate;
    /**
     * \brief EPA*SE relevant data
     */
    unsigned int v;
    /**
     * \brief EPA*SE relevant data
     */
    unsigned int g;
    /**
     * \brief EPA*SE relevant data
     */
    unsigned int gp;
    /**
     * \brief EPA*SE relevant data
     */
    short unsigned int iterationclosed;
    /**
     * \brief EPA*SE relevant data
     */
    short unsigned int callnumberaccessed;

#if DEBUG
    /**
     * \brief EPA*SE relevant data
     */
    short unsigned int numofexpands;
#endif

    /**
     * \brief EPA*SE relevant data
     */
    short unsigned int mask;
    /**
     * \brief EPA*SE relevant data
     */
    std::multimap<int,EPASESearchStateData*>::const_iterator iter;
    /**
     * \brief best predecessor and the action from it, used only in forward searches
     */
    CMDPSTATE *bestpredstate;
    /**
     * \brief the next state if executing best action
     */
    CMDPSTATE *bestnextstate;
    unsigned int costtobestnextstate;
    unsigned int h;

public:
    EPASESearchStateData() { }
    ~EPASESearchStateData() { }
} EPASEState;

class myHeap { 
  private:
    std::multimap<int,EPASEState*> open, be;
    bool* done;
    DiscreteSpaceInformation* env;
    double eps;
    bool usable;
  public:
    myHeap(bool* done_flag, DiscreteSpaceInformation* e);
    void setEps(double e);
    int bound(const EPASEState* const s) const;
    void clear();
    bool empty() const;
    void initIter(EPASEState* const s) const;
    void insert(EPASEState* const s);
    void removeBE(EPASEState* const s);
    int min_value();
    void analyze();
    EPASEState* remove(std::unique_lock<std::mutex>* lock, int* fval, int thread_id);

    class DualIterator {
      public:
        std::multimap<int,EPASEState*>::const_iterator a_iter, b_iter;
        DualIterator(const myHeap& q) : a_iter(q.open.cbegin()), b_iter(q.be.cbegin()) {}
        EPASEState* next() {
          if (a_iter->first < b_iter->first) {
            return (a_iter++)->second;
          } else {
            return (b_iter++)->second;
          }
        }
    };
};

/**
 * \brief the statespace of EPA*SE
 */
typedef struct EPASESearchStateSpace
{
    double eps;
    double eps_satisfied;
    myHeap* heap;
    CList* inconslist;
    short unsigned int searchiteration;
    short unsigned int callnumber;
    CMDPSTATE* searchgoalstate;
    CMDPSTATE* searchstartstate;

    CMDP searchMDP;

    bool bReevaluatefvals;
    bool bReinitializeSearchStateSpace;
    bool bNewSearchIteration;
} EPASESearchStateSpace_t;

/**
 * \brief EPASE* planner
 */
class EPASEPlanner : public SBPLPlanner
{
  public:
    /**
     * \brief replan a path within the allocated time, return the solution in the vector
     */
    virtual int replan(double allocated_time_secs, std::vector<int>* solution_stateIDs_V);

    /**
     * \brief replan a path within the allocated time, return the solution in the vector, also returns solution cost
     */
    virtual int replan(double allocated_time_sec, std::vector<int>* solution_stateIDs_V, int* solcost);

    /**
     * \brief works same as replan function with time and solution states, but
     *        it let's you fill out all the parameters for the search
     */
    virtual int replan(std::vector<int>* solution_stateIDs_V, ReplanParams params);

    /**
     * \brief works same as replan function with time, solution states, and
     *        cost, but it let's you fill out all the parameters for the search
     */
    virtual int replan(std::vector<int>* solution_stateIDs_V, ReplanParams params, int* solcost);

    /**
     * \brief set the goal state
     */
    virtual int set_goal(int goal_stateID);

    /**
     * \brief set the start state
     */
    virtual int set_start(int start_stateID);

    /**
     * \brief inform the search about the new edge costs
     */
    virtual void costs_changed(StateChangeQuery const & stateChange);

    /**
     * \brief inform the search about the new edge costs -
     * \note since EPASE* is non-incremental, it is sufficient (and more
     *       efficient) to just inform EPASE* of the fact that some costs changed
     */
    virtual void costs_changed();

    /**
     * \brief set a flag to get rid of the previous search efforts, and
     *        re-initialize the search, when the next replan is called
     */
    virtual int force_planning_from_scratch();

    /**
     * \brief Gets rid of the previous search efforts, release the memory and re-initialize the search.
     */
    virtual int force_planning_from_scratch_and_free_memory();

    /**
     * \brief you can either search forwards or backwards
     */
    virtual int set_search_mode(bool bSearchUntilFirstSolution);

    /**
     * \brief returns the suboptimality bound on the currently found solution
     */
    virtual double get_solution_eps() const { return pSearchStateSpace_->eps_satisfied; }

    /**
     * \brief returns the number of states expanded so far
     */
    virtual int get_n_expands() const { return searchexpands; }

    /**
     * \brief sets the value of the initial epsilon (suboptimality bound) used
     */
    virtual void set_initialsolution_eps(double initialsolution_eps) { finitial_eps = initialsolution_eps; }

    /**
     * \brief sets the value to decrease from eps at each iteration
     */
    virtual void set_eps_step(double eps) { dec_eps = eps; }

    /**
     * \brief prints out the search path into a file
     */
    virtual void print_searchpath(FILE* fOut);

    /**
     * \brief constructor
     */
    EPASEPlanner(DiscreteSpaceInformation* environment, bool bforwardsearch);

    /**
     * \brief destructor
     */
    ~EPASEPlanner();

    virtual void astarThread();

    /**
     * \brief returns the initial epsilon
     */
    virtual double get_initial_eps() { return finitial_eps; }

    /**
     * \brief returns the time taken to find the first solution
     */
    virtual double get_initial_eps_planning_time() { return finitial_eps_planning_time; }

    /**
     * \brief returns the time taken to get the final solution
     */
    virtual double get_final_eps_planning_time() { return final_eps_planning_time; }

    /**
     * \brief Return the number of expands to find the first solution or -1 if no solution has been found.
     */
    virtual int get_n_expands_init_solution() { return num_of_expands_initial_solution; }

    /**
     * \brief returns the final epsilon achieved during the search
     */
    virtual double get_final_epsilon() { return final_eps; }

    /**
     * \brief fills out a vector of stats from the search
     */
    virtual void get_search_stats(std::vector<PlannerStats>* s);

  protected:
    //member variables
    double finitial_eps, finitial_eps_planning_time, final_eps_planning_time, final_eps, dec_eps, final_epsilon;
    double repair_time;
    bool use_repair_time;

    std::vector<PlannerStats> stats;

    int num_of_expands_initial_solution;

    MDPConfig* MDPCfg_;

    bool bforwardsearch; //if true, then search proceeds forward, otherwise backward

    bool bsearchuntilfirstsolution; //if true, then search until first solution only (see planner.h for search modes)

    EPASESearchStateSpace_t* pSearchStateSpace_;

    unsigned int searchexpands;
    int MaxMemoryCounter;
    clock_t TimeStarted;
    FILE *fDeb;

    // parallel stuff
    std::mutex the_mutex;
    std::condition_variable worker_cond;
    std::condition_variable main_cond;
    bool planner_ok;
    bool iteration_done;
    int sleeping_threads;
    int thread_ids;
    std::vector<std::thread> threads;
    int improve_path_result;
    int bad_cnt;

    //member functions
    virtual void Initialize_searchinfo(CMDPSTATE* state, EPASESearchStateSpace_t* pSearchStateSpace);

    virtual CMDPSTATE* CreateState(int stateID, EPASESearchStateSpace_t* pSearchStateSpace);

    virtual CMDPSTATE* GetState(int stateID, EPASESearchStateSpace_t* pSearchStateSpace);

    virtual int ComputeHeuristic(CMDPSTATE* MDPstate, EPASESearchStateSpace_t* pSearchStateSpace);

    //initialization of a state
    virtual void InitializeSearchStateInfo(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace);

    //re-initialization of a state
    virtual void ReInitializeSearchStateInfo(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace);

    virtual void DeleteSearchStateData(EPASEState* state);

    //used for backward search
    virtual void UpdatePreds(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace);

    //used for forward search
    virtual void UpdateSuccs(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace);

    virtual int GetGVal(int StateID, EPASESearchStateSpace_t* pSearchStateSpace);

    //returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
    virtual int ImprovePath(EPASESearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);

    virtual void BuildNewOPENList(EPASESearchStateSpace_t* pSearchStateSpace);

    virtual void Reevaluatefvals(EPASESearchStateSpace_t* pSearchStateSpace);
    virtual void Reevaluatehvals(EPASESearchStateSpace_t* pSearchStateSpace);

    //creates (allocates memory) search state space
    //does not initialize search statespace
    virtual int CreateSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace);

    //deallocates memory used by SearchStateSpace
    virtual void DeleteSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace);

    //debugging
    virtual void PrintSearchState(EPASEState* state, FILE* fOut);

    //reset properly search state space
    //needs to be done before deleting states
    virtual int ResetSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace);

    //initialization before each search
    virtual void ReInitializeSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace);

    //very first initialization
    virtual int InitializeSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace);

    virtual int SetSearchGoalState(int SearchGoalStateID, EPASESearchStateSpace_t* pSearchStateSpace);

    virtual int SetSearchStartState(int SearchStartStateID, EPASESearchStateSpace_t* pSearchStateSpace);

    //reconstruct path functions are only relevant for forward search
    virtual int ReconstructPath(EPASESearchStateSpace_t* pSearchStateSpace);

    virtual void PrintSearchPath(EPASESearchStateSpace_t* pSearchStateSpace, FILE* fOut);

    virtual int getHeurValue(EPASESearchStateSpace_t* pSearchStateSpace, int StateID);

    //get path
    virtual std::vector<int> GetSearchPath(EPASESearchStateSpace_t* pSearchStateSpace, int& solcost);

    virtual bool Search(EPASESearchStateSpace_t* pSearchStateSpace, std::vector<int>& pathIds, int & PathCost,
                        bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs);

    virtual bool outOfTime();
};

#endif

