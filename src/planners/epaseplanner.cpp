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

#include <cmath>
#include <sbpl/discrete_space_information/environment.h>
#include <sbpl/planners/epaseplanner.h>
#include <sbpl/utils/heap.h>
#include <sbpl/utils/key.h>
#include <sbpl/utils/list.h>

using namespace std;

//-------------------------------------------------------------------------------------------
//Auxiliary functions
typedef chrono::high_resolution_clock Clock;

double wA_eps, pA_eps, leeway; // TODO: initialize
Clock::time_point timeStart;

double GetTimeSince(const Clock::time_point& ts) {
  Clock::time_point t = Clock::now();
  double td = chrono::duration<double, chrono::seconds::period>(t-ts).count();
  return td;
}

inline int bound_back(const EPASEState* const sp, const EPASEState* const s) {
  if (wA_eps <= pA_eps)
    return s->g + sp->iter->first - s->iter->first + leeway;
  else
    return (s->g + sp->iter->first - s->iter->first)*pA_eps/wA_eps + leeway;
}


//-------------------------------------------------------------------------------------------
//OPEN list
myHeap::myHeap(bool* done_flag, DiscreteSpaceInformation* e){
  done = done_flag;
  env = e;
}

void myHeap::setEps(double e){
  eps = e;
}

int myHeap::bound(const EPASEState* const s) const {
  myHeap::DualIterator q_it(*this);
  int g_front = s->gp;
  EPASEState* sp = q_it.next();
  int g_back = bound_back(sp, s);
  while (g_back < s->g && s->g <= g_front) {
    int next_cost = sp->gp + pA_eps*env->fastHeuristic(sp->id,s->id);
    g_front = min(g_front, next_cost);
    sp = q_it.next();
    g_back = bound_back(sp, s);
  }
  return min(g_front, g_back);
}

void myHeap::clear(){
  open.clear(); be.clear(); /*closed.clear();
  null_state.mask = null_state.x = null_state.y = null_state.gp = 0;
  null_state.g = INFINITE;
  null_state.iter =
    open.insert(pair<Cost,State*>(f(&null_state),&null_state));
      be.insert(pair<Cost,State*>(f(&null_state),&null_state));*/
}

bool myHeap::empty() const {
  return open.empty();
}

void myHeap::initIter(EPASEState* const s) const {
  s->iter = open.cend();
}

void myHeap::insert(EPASEState* const s){
  if (s->iter != open.cend())
    open.erase(s->iter);
  s->iter = open.insert(pair<int,EPASEState*>(s->g + wA_eps*ComputeHeuristic(s->MDPstate, a), s));
}

void myHeap::removeBE(EPASEState* const s){
  be.erase(s->iter);
  initIter(s);
}

int myHeap::min_value(){
  return open.begin()->first;
}

void myHeap::analyze(){
  printf("queue={\n");
  int cnt = 0;
  for(multimap<int,EPASEState*>::const_iterator it=open.cbegin(); it!=open.cend(); it++){
    if(cnt > 10){
      printf("...\n");
      break;
    }
    printf("       f=%d, g=%d, h=%d: \n",it->first,it->second->g,env->fastHeuristic(open.cbegin()->second->id,it->second->id));
    //env->PrintState(it->second->id,true);
    cnt++;
  }
  printf("}\n");
  int ties = -1;
  int first_f = open.cbegin()->first;
  for(multimap<int,EPASEState*>::const_iterator it=open.cbegin(); it!=open.cend(); it++){
    if(it->first == first_f)
      ties++;
  }
  vector<int> validV;
  cnt = 0;
  for(multimap<int,EPASEState*>::const_iterator it=open.cbegin(); it!=open.cend(); it++){
    if(cnt > 10){
      printf("...\n");
      break;
    }
    bool valid = true;
    int cnt2 = 0;
    for(multimap<int,EPASEState*>::const_iterator it2=open.cbegin(); it2!=it; it2++){
      if(int(it->second->g) - int(it2->second->g) > eps*env->fastHeuristic(it->second->id,it2->second->id)){
        printf("invalid (%d,%d): %d - %d > %f*%d\n",
                cnt,cnt2,
                it->second->g, it2->second->g, eps, 
                env->fastHeuristic(it->second->id,it2->second->id));
        valid = false;
        break;
      }
      cnt2++;
    }
    if(valid)
      validV.push_back(cnt);
    cnt++;
  }
  printf("parallel states={");
  for(unsigned int i=0; i<validV.size(); i++)
    printf("%d,",validV[i]);
  printf("}\n");
  printf("%d parallel states with %d f-val ties\n",int(validV.size()),ties);
  //cin.get();
}

EPASEState* myHeap::remove(int& gb, int& fval){
  //analyze();
  
  //loop over open list to find a State we can expand
  if (usable) {
    for (multimap<int,pARAState*>::const_iterator it = open.cbegin(); it != open.cend(); ++it) {
      State* s = it->second;
      gb = bound(s);
      if (s->g <= gb) {
        //printf("got a State!");
        fval = it->first;
        open.erase(s->iter);
        s->iter = be.insert(pair<int,pARAState*>(f(s),s));
        /*#ifdef ANYTIME
          closed.push_back(s);
        #endif
        */
        return s;
      }
    }
  }
  usable = false;
  return NULL;
}

//-------------------------------------------------------------------------------------------

EPASEPlanner::EPASEPlanner(DiscreteSpaceInformation* environment, bool bSearchForward) :
  heap(&iteration_done, environment), params(0.0) {
  bforwardsearch = bSearchForward;
  environment_ = environment;
  replan_number = 0;
  //reconstructTime = 0.0;

  bsearchuntilfirstsolution = false;
  finitial_eps = EPASE_DEFAULT_INITIAL_EPS;
  final_epsilon = EPASE_FINAL_EPS;
  dec_eps = EPASE_DECREASE_EPS;
  use_repair_time = false;
  repair_time = INFINITECOST;
  searchexpands = 0;
  MaxMemoryCounter = 0;

#ifndef ROS
  const char* debug = "debug.txt";
#endif
  fDeb = SBPL_FOPEN(debug, "w");
  if (fDeb == NULL) {
    SBPL_ERROR("ERROR: could not open planner debug file\n");
    throw new SBPL_Exception();
  }

  pSearchStateSpace_ = new EPASESearchStateSpace_t;

  //create the EPASE planner
  if (CreateSearchStateSpace(pSearchStateSpace_) != 1) {
    SBPL_ERROR("ERROR: failed to create statespace\n");
    return;
  }

  //set the start and goal states
  if (InitializeSearchStateSpace(pSearchStateSpace_) != 1) {
    SBPL_ERROR("ERROR: failed to create statespace\n");
    return;
  }
  finitial_eps_planning_time = -1.0;
  final_eps_planning_time = -1.0;
  num_of_expands_initial_solution = 0;
  final_eps = -1.0;

  fout = fopen("stats.txt", "w");
  // TODO: push start state

  goal_state_id = -1;
  start_state_id = -1;

  unique_lock<mutex> lock(mutex);
  //vector<boost::thread*> threads;
  //vector<pARAState*> being_expanded;
  planner_ok = true;
  thread_ids = 0;
  for(int i=0; i<EPASE_NUM_THREADS; i++){
    threads.push_back(thread(&EPASEPlanner::astarThread));
  }
  main_cond.wait(lock);
}

EPASEPlanner::~EPASEPlanner() {
  fclose(fout);
  unique_lock<mutex> lock(the_mutex);
  planner_ok = false;
  worker_cond.notify_all();
  lock.unlock();
  for(int i=0; i<EPASE_NUM_THREADS; i++)
    threads[i].join();

  if (pSearchStateSpace_ != NULL) {
    //delete the statespace
    DeleteSearchStateSpace( pSearchStateSpace_);
    delete pSearchStateSpace_;
  }
  SBPL_FCLOSE( fDeb);
}

void EPASEPlanner::astarThread() {
  unique_lock<mutex> thread_lock(mutex);
  int thread_id = thread_ids;
  thread_ids++;
  while (true) { // TODO find wakeup condition and then main_cond.notify_one();
    //the mutex is locked
      
    worker_cond.wait(thread_lock);
    thread_lock.unlock();
    if (!planner_ok)
      break;
    
    int ret = ImprovePath(thread_id);
    if (ret >= 0){
      printf("thread %d: setting the improve_path_result to %d\n",thread_id,ret);
      improve_path_result = ret;
    }
    thread_lock.lock();
  }
}

void EPASEPlanner::Initialize_searchinfo(CMDPSTATE* state, EPASESearchStateSpace_t* pSearchStateSpace)
{
    EPASEState* searchstateinfo = (EPASEState*)state->PlannerSpecificData;

    searchstateinfo->MDPstate = state;
    InitializeSearchStateInfo(searchstateinfo, pSearchStateSpace);
}

CMDPSTATE* EPASEPlanner::CreateState(int stateID, EPASESearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* state = NULL;

#if DEBUG
    if (environment_->StateID2IndexMapping[stateID][EPASEMDP_STATEID2IND] != -1) {
        SBPL_ERROR("ERROR in CreateState: state already created\n");
        throw new SBPL_Exception();
    }
#endif

    //adds to the tail a state
    state = pSearchStateSpace->searchMDP.AddState(stateID);

    //remember the index of the state
    environment_->StateID2IndexMapping[stateID][EPASEMDP_STATEID2IND] =
            pSearchStateSpace->searchMDP.StateArray.size() - 1;

#if DEBUG
    if(state !=
       pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][EPASEMDP_STATEID2IND]])
    {
        SBPL_ERROR("ERROR in CreateState: invalid state index\n");
        throw new SBPL_Exception();
    }
#endif

    //create search specific info
    state->PlannerSpecificData = (EPASEState*)malloc(sizeof(EPASEState));
    Initialize_searchinfo(state, pSearchStateSpace);
    MaxMemoryCounter += sizeof(EPASEState);

    return state;
}

CMDPSTATE* EPASEPlanner::GetState(int stateID, EPASESearchStateSpace_t* pSearchStateSpace)
{
    if (stateID >= (int)environment_->StateID2IndexMapping.size()) {
        SBPL_ERROR("ERROR int GetState: stateID %d is invalid\n", stateID);
        throw new SBPL_Exception();
    }

    if (environment_->StateID2IndexMapping[stateID][EPASEMDP_STATEID2IND] == -1)
        return CreateState(stateID, pSearchStateSpace);
    else
        return pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][EPASEMDP_STATEID2IND]];
}

//-----------------------------------------------------------------------------------------------------

int EPASEPlanner::ComputeHeuristic(CMDPSTATE* MDPstate, EPASESearchStateSpace_t* pSearchStateSpace)
{
    //compute heuristic for search

    if (bforwardsearch) {

#if MEM_CHECK == 1
        //int WasEn = DisableMemCheck();
#endif

        //forward search: heur = distance from state to searchgoal which is Goal EPASEState
        int retv = environment_->GetGoalHeuristic(MDPstate->StateID);

#if MEM_CHECK == 1
        //if (WasEn)
        //	EnableMemCheck();
#endif

        return retv;

    }
    else {
        //backward search: heur = distance from searchgoal to state
        return environment_->GetStartHeuristic(MDPstate->StateID);
    }
}

// initialization of a state
void EPASEPlanner::InitializeSearchStateInfo(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace)
{
    state->g = INFINITECOST;
    state->v = INFINITECOST;
    state->iterationclosed = 0;
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->bestnextstate = NULL;
    state->costtobestnextstate = INFINITECOST;
    state->heapindex = 0;
    state->listelem[EPASE_INCONS_LIST_ID] = 0;
#if DEBUG
    state->numofexpands = 0;
#endif

    state->bestpredstate = NULL;

    //compute heuristics
#if USE_HEUR
    if(pSearchStateSpace->searchgoalstate != NULL)
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    else
        state->h = 0;
#else
        state->h = 0;
#endif
}

//re-initialization of a state
void EPASEPlanner::ReInitializeSearchStateInfo(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace)
{
    state->g = INFINITECOST;
    state->v = INFINITECOST;
    state->iterationclosed = 0;
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->bestnextstate = NULL;
    state->costtobestnextstate = INFINITECOST;
    state->heapindex = 0;
    state->listelem[EPASE_INCONS_LIST_ID] = 0;
#if DEBUG
    state->numofexpands = 0;
#endif

    state->bestpredstate = NULL;

    //compute heuristics
#if USE_HEUR
    if(pSearchStateSpace->searchgoalstate != NULL)
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    else
        state->h = 0;
#else
        state->h = 0;
#endif
}

void EPASEPlanner::DeleteSearchStateData(EPASEState* state)
{
    //no memory was allocated
    MaxMemoryCounter = 0;
    return;
}

//used for backward search
void EPASEPlanner::UpdatePreds(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace)
{
    vector<int> PredIDV;
    vector<int> CostV;
    CKey key;
    EPASEState *predstate;

    environment_->GetPreds(state->MDPstate->StateID, &PredIDV, &CostV);

    //iterate through predecessors of s
    for (int pind = 0; pind < (int)PredIDV.size(); pind++) {
        CMDPSTATE* PredMDPState = GetState(PredIDV[pind], pSearchStateSpace);
        predstate = (EPASEState*)(PredMDPState->PlannerSpecificData);
        if (predstate->callnumberaccessed != pSearchStateSpace->callnumber) {
            ReInitializeSearchStateInfo(predstate, pSearchStateSpace);
        }

        //see if we can improve the value of predstate
        if (predstate->g > state->v + CostV[pind]) {
            predstate->g = state->v + CostV[pind];
            predstate->bestnextstate = state->MDPstate;
            predstate->costtobestnextstate = CostV[pind];

            //re-insert into heap if not closed yet
            if (predstate->iterationclosed != pSearchStateSpace->searchiteration) {
                key.key[0] = predstate->g + (int)(pSearchStateSpace->eps * predstate->h);
                //key.key[1] = predstate->h;
                if (predstate->heapindex != 0)
                    pSearchStateSpace->heap->updateheap(predstate, key);
                else
                    pSearchStateSpace->heap->insertheap(predstate, key);
            }
            //take care of incons list
            else if (predstate->listelem[EPASE_INCONS_LIST_ID] == NULL) {
                pSearchStateSpace->inconslist->insert(predstate, EPASE_INCONS_LIST_ID);
            }
        }
    } //for predecessors
}

//used for forward search
void EPASEPlanner::UpdateSuccs(EPASEState* state, EPASESearchStateSpace_t* pSearchStateSpace)
{
    vector<int> SuccIDV;
    vector<int> CostV;
    CKey key;
    EPASEState *succstate;

    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);

    //iterate through predecessors of s
    for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
        CMDPSTATE* SuccMDPState = GetState(SuccIDV[sind], pSearchStateSpace);
        int cost = CostV[sind];

        succstate = (EPASEState*)(SuccMDPState->PlannerSpecificData);
        if (succstate->callnumberaccessed != pSearchStateSpace->callnumber) {
            ReInitializeSearchStateInfo(succstate, pSearchStateSpace);
        }

        //see if we can improve the value of succstate
        //taking into account the cost of action
        if (succstate->g > state->v + cost) {
            succstate->g = state->v + cost;
            succstate->bestpredstate = state->MDPstate;

            //re-insert into heap if not closed yet
            if (succstate->iterationclosed != pSearchStateSpace->searchiteration) {

                key.key[0] = succstate->g + (int)(pSearchStateSpace->eps * succstate->h);

                //key.key[1] = succstate->h;

                if (succstate->heapindex != 0)
                    pSearchStateSpace->heap->updateheap(succstate, key);
                else
                    pSearchStateSpace->heap->insertheap(succstate, key);
            }
            //take care of incons list
            else if (succstate->listelem[EPASE_INCONS_LIST_ID] == NULL) {
                pSearchStateSpace->inconslist->insert(succstate, EPASE_INCONS_LIST_ID);
            }
        } //check for cost improvement
    } //for actions
}

//TODO-debugmax - add obsthresh and other thresholds to other environments in 3dkin
int EPASEPlanner::GetGVal(int StateID, EPASESearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* cmdp_state = GetState(StateID, pSearchStateSpace);
    EPASEState* state = (EPASEState*)cmdp_state->PlannerSpecificData;
    return state->g;
}

//returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
int EPASEPlanner::ImprovePath(EPASESearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs)
{
    int expands, fval, gb;
    EPASEState *state, *searchgoalstate;
    CKey key, minkey;
    CKey goalkey;

    expands = 0;

    if (pSearchStateSpace->searchgoalstate == NULL) {
        SBPL_ERROR("ERROR searching: no goal state is set\n");
        throw new SBPL_Exception();
    }

    //goal state
    searchgoalstate = (EPASEState*)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
    if (searchgoalstate->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo(searchgoalstate, pSearchStateSpace);
    }

    //set goal key
    goalkey.key[0] = searchgoalstate->g;
    //goalkey.key[1] = searchgoalstate->h;

    //expand states until done
    minkey = pSearchStateSpace->heap.remove(gb, fval);
    CKey oldkey = minkey;
    while (!pSearchStateSpace->heap.empty() && minkey.key[0] < INFINITECOST && goalkey > minkey &&
           GetTimeSince(timeStart) < MaxNumofSecs &&
               (pSearchStateSpace->eps_satisfied == INFINITECOST ||
               GetTimeSince(timeStart) < repair_time ))
    {
		//get the state
		state = (EPASEState*)pSearchStateSpace->heap.remove(gb, fval);
		/* TODO: look at this old epase code:
		pARAState* state = heap.remove(&lock, &fval,thread_id);
		if(state==NULL){
			printf("thread %d: heap remove returned NULL\n",thread_id);
			break;
		}
    */

#if DEBUG
        SBPL_FPRINTF(fDeb, "expanding state(%d): h=%d g=%u key=%u v=%u iterc=%d callnuma=%d expands=%d (g(goal)=%u)\n",
                     state->MDPstate->StateID, state->h, state->g, state->g+(int)(pSearchStateSpace->eps*state->h),
                     state->v, state->iterationclosed, state->callnumberaccessed, state->numofexpands,
                     searchgoalstate->g);
        SBPL_FPRINTF(fDeb, "expanding: ");
        PrintSearchState(state, fDeb);
        if (state->listelem[EPASE_INCONS_LIST_ID] != NULL) {
            SBPL_FPRINTF(fDeb, "ERROR: expanding a state from inconslist\n");
            SBPL_ERROR("ERROR: expanding a state from inconslist\n");
            throw new SBPL_Exception();
        }
        //SBPL_FFLUSH(fDeb);
#endif

#if DEBUG
        if (minkey.key[0] < oldkey.key[0] && fabs(this->finitial_eps - 1.0) < ERR_EPS) {
            //SBPL_PRINTF("WARN in search: the sequence of keys decreases\n");
            //throw new SBPL_Exception();
        }
        oldkey = minkey;
#endif

		if(goal_state->g <= (unsigned int)fval){
			printf("thread %d: goal g-val is less than min fval\n",thread_id);
			iteration_done = true;
			return 1;
			break;
		}

        if (state->v == state->g) {
            SBPL_ERROR("ERROR: consistent state is being expanded\n");
#if DEBUG
            SBPL_FPRINTF(fDeb, "ERROR: consistent state is being expanded\n");
            throw new SBPL_Exception();
#endif
        }

        pSearchStateSpace->heap.removeBE(state);

        //recompute state value
        state->v = state->g;
        state->iterationclosed = pSearchStateSpace->searchiteration;

        //new expand
        expands++;
#if DEBUG
        state->numofexpands++;
#endif

    int cnt = 1;
    //cnt += q.be.size();
    //printf("thread %d: expanding with %d other threads\n",thread_id,cnt);
    //environment_->PrintState(state->id,true);
    if(thread_id==0){
      rate += cnt;
      rate_cnt++;
      if(rate_cnt==100){
        printf("avg parallel threads = %f\n", rate/rate_cnt);
        rate = 0;
        rate_cnt = 0;
      }
    }
    if(cnt==0)
      bad_cnt++;
    else
      bad_cnt=0;
    //if(bad_cnt>3)
      //heap.analyze();
    being_expanded[thread_id] = state;
    //being_expanded_fval[thread_id] = fval;
    lock.unlock();

        if (bforwardsearch)
            UpdateSuccs(state, pSearchStateSpace);
        else
            UpdatePreds(state, pSearchStateSpace);

        //recompute minkey
        minkey = pSearchStateSpace->heap->getminkeyheap();

        //recompute goalkey if necessary
        if (goalkey.key[0] != (int)searchgoalstate->g) {
            //SBPL_PRINTF("re-computing goal key\n");
            //recompute the goal key (heuristics should be zero)
            goalkey.key[0] = searchgoalstate->g;
            //goalkey.key[1] = searchgoalstate->h;
        }

        if (expands % 100000 == 0 && expands > 0) {
            SBPL_PRINTF("expands so far=%u\n", expands);
        }
    }

    int retv = 1;
    if (searchgoalstate->g == INFINITECOST && pSearchStateSpace->heap->emptyheap()) {
        SBPL_PRINTF("solution does not exist: search exited because heap is empty\n");
        retv = 0;
    }
    else if (!pSearchStateSpace->heap->emptyheap() && goalkey > minkey) {
        SBPL_PRINTF("search exited because it ran out of time\n");
        retv = 2;
    }
    else if (searchgoalstate->g == INFINITECOST && !pSearchStateSpace->heap->emptyheap()) {
        SBPL_PRINTF("solution does not exist: search exited because all candidates for expansion have "
                    "infinite heuristics\n");
        retv = 0;
    }
    else {
        SBPL_PRINTF("search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
        retv = 1;
    }

    //SBPL_FPRINTF(fDeb, "expanded=%d\n", expands);

    searchexpands += expands;

    return retv;
}

void EPASEPlanner::BuildNewOPENList(EPASESearchStateSpace_t* pSearchStateSpace)
{
    EPASEState *state;
    CKey key;
    myHeap* pheap = pSearchStateSpace->heap;
    CList* pinconslist = pSearchStateSpace->inconslist;

    //move incons into open
    while (pinconslist->firstelement != NULL) {
        state = (EPASEState*)pinconslist->firstelement->liststate;

        //compute f-value
        key.key[0] = state->g + (int)(pSearchStateSpace->eps * state->h);
        //key.key[1] = state->h;

        //insert into OPEN
        pheap->insertheap(state, key);
        //remove from INCONS
        pinconslist->remove(state, EPASE_INCONS_LIST_ID);
    }
}

void EPASEPlanner::Reevaluatefvals(EPASESearchStateSpace_t* pSearchStateSpace)
{
    CKey key;
    int i;
    myHeap* pheap = pSearchStateSpace->heap;

    //recompute priorities for states in OPEN and reorder it
    for (i = 1; i <= pheap->currentsize; ++i) {
        EPASEState* state = (EPASEState*)pheap->heap[i].heapstate;
        pheap->heap[i].key.key[0] = state->g + (int)(pSearchStateSpace->eps * state->h);
        //pheap->heap[i].key.key[1] = state->h;
    }
    pheap->makeheap();

    pSearchStateSpace->bReevaluatefvals = false;
}

void EPASEPlanner::Reevaluatehvals(EPASESearchStateSpace_t* pSearchStateSpace)
{
    for(int i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++)
    {
        CMDPSTATE* MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
        EPASEState* state = (EPASEState*)MDPstate->PlannerSpecificData;
        state->h = ComputeHeuristic(MDPstate, pSearchStateSpace);
    }
}

//creates (allocates memory) search state space
//does not initialize search statespace
int EPASEPlanner::CreateSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace)
{
    //create a heap
    pSearchStateSpace->heap = new myHeap;
    pSearchStateSpace->inconslist = new CList;
    MaxMemoryCounter += sizeof(myHeap);
    MaxMemoryCounter += sizeof(CList);

    pSearchStateSpace->searchgoalstate = NULL;
    pSearchStateSpace->searchstartstate = NULL;

    searchexpands = 0;
    num_of_expands_initial_solution = -1;

    pSearchStateSpace->bReinitializeSearchStateSpace = false;

    return 1;
}

//deallocates memory used by SearchStateSpace
void EPASEPlanner::DeleteSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->heap != NULL) {
        pSearchStateSpace->heap->makeemptyheap();
        delete pSearchStateSpace->heap;
        pSearchStateSpace->heap = NULL;
    }

    if (pSearchStateSpace->inconslist != NULL) {
        pSearchStateSpace->inconslist->makeemptylist(EPASE_INCONS_LIST_ID);
        delete pSearchStateSpace->inconslist;
        pSearchStateSpace->inconslist = NULL;
    }

    //delete the states themselves
    int iend = (int)pSearchStateSpace->searchMDP.StateArray.size();
    for (int i = 0; i < iend; i++) {
        CMDPSTATE* state = pSearchStateSpace->searchMDP.StateArray[i];
        if (state != NULL && state->PlannerSpecificData != NULL) {
            DeleteSearchStateData((EPASEState*)state->PlannerSpecificData);
            free((EPASEState*)state->PlannerSpecificData);
            state->PlannerSpecificData = NULL;
        }
    }
    pSearchStateSpace->searchMDP.Delete();
}

//reset properly search state space
//needs to be done before deleting states
int EPASEPlanner::ResetSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace)
{
    pSearchStateSpace->heap->makeemptyheap();
    pSearchStateSpace->inconslist->makeemptylist(EPASE_INCONS_LIST_ID);

    return 1;
}

//initialization before each search
void EPASEPlanner::ReInitializeSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace)
{
    CKey key;

    //increase callnumber
    pSearchStateSpace->callnumber++;

    //reset iteration
    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->bNewSearchIteration = true;

#if DEBUG
    SBPL_FPRINTF(fDeb, "reinitializing search state-space (new call number=%d search iter=%d)\n",
        pSearchStateSpace->callnumber,pSearchStateSpace->searchiteration );
#endif

    pSearchStateSpace->heap->makeemptyheap();
    pSearchStateSpace->inconslist->makeemptylist(EPASE_INCONS_LIST_ID);

    //reset 
    pSearchStateSpace->eps = this->finitial_eps;
    pSearchStateSpace->eps_satisfied = INFINITECOST;

    //initialize start state
    EPASEState* startstateinfo = (EPASEState*)(pSearchStateSpace->searchstartstate->PlannerSpecificData);
    if (startstateinfo->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo(startstateinfo, pSearchStateSpace);
    }
    startstateinfo->g = 0;
    
    //initialize goal state
    EPASEState* searchgoalstate = (EPASEState*)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
    if (searchgoalstate->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo(searchgoalstate, pSearchStateSpace);
    }

    //insert start state into the heap
    key.key[0] = (long int)(pSearchStateSpace->eps * startstateinfo->h);
    //key.key[1] = startstateinfo->h;
    pSearchStateSpace->heap->insertheap(startstateinfo, key);

    pSearchStateSpace->bReinitializeSearchStateSpace = false;
    pSearchStateSpace->bReevaluatefvals = false;
}

//very first initialization
int EPASEPlanner::InitializeSearchStateSpace(EPASESearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->heap->currentsize != 0 || pSearchStateSpace->inconslist->currentsize != 0) {
        SBPL_ERROR("ERROR in InitializeSearchStateSpace: heap or list is not empty\n");
        throw new SBPL_Exception();
    }

    pSearchStateSpace->eps = this->finitial_eps;
    pSearchStateSpace->eps_satisfied = INFINITECOST;
    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->bNewSearchIteration = true;
    pSearchStateSpace->callnumber = 0;
    pSearchStateSpace->bReevaluatefvals = false;

    //create and set the search start state
    pSearchStateSpace->searchgoalstate = NULL;
    //pSearchStateSpace->searchstartstate = GetState(SearchStartStateID, pSearchStateSpace);
    pSearchStateSpace->searchstartstate = NULL;

    pSearchStateSpace->bReinitializeSearchStateSpace = true;

    return 1;
}

int EPASEPlanner::SetSearchGoalState(int SearchGoalStateID, EPASESearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->searchgoalstate == NULL ||
        pSearchStateSpace->searchgoalstate->StateID != SearchGoalStateID)
    {
        pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);

        //should be new search iteration
        pSearchStateSpace->eps_satisfied = INFINITECOST;
        pSearchStateSpace->bNewSearchIteration = true;
        pSearchStateSpace_->eps = this->finitial_eps;

#if USE_HEUR
        //recompute heuristic for the heap if heuristics is used
        pSearchStateSpace->bReevaluatefvals = true;
#endif
    }

    return 1;
}

int EPASEPlanner::SetSearchStartState(int SearchStartStateID, EPASESearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* MDPstate = GetState(SearchStartStateID, pSearchStateSpace);

    if (MDPstate != pSearchStateSpace->searchstartstate) {
        pSearchStateSpace->searchstartstate = MDPstate;
        pSearchStateSpace->bReinitializeSearchStateSpace = true;
    }

    return 1;
}

int EPASEPlanner::ReconstructPath(EPASESearchStateSpace_t* pSearchStateSpace)
{
    if (bforwardsearch) //nothing to do, if search is backward
    {
        CMDPSTATE* MDPstate = pSearchStateSpace->searchgoalstate;
        CMDPSTATE* PredMDPstate;
        EPASEState *predstateinfo, *stateinfo;

#if DEBUG
        SBPL_FPRINTF(fDeb, "reconstructing a path:\n");
#endif

        while (MDPstate != pSearchStateSpace->searchstartstate) {
            stateinfo = (EPASEState*)MDPstate->PlannerSpecificData;

#if DEBUG
            PrintSearchState(stateinfo, fDeb);
#endif
            if (stateinfo->g == INFINITECOST) {
                //SBPL_ERROR("ERROR in ReconstructPath: g of the state on the path is INFINITE\n");
                //throw new SBPL_Exception();
                return -1;
            }

            if (stateinfo->bestpredstate == NULL) {
                SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
                throw new SBPL_Exception();
            }

            //get the parent state
            PredMDPstate = stateinfo->bestpredstate;
            predstateinfo = (EPASEState*)PredMDPstate->PlannerSpecificData;

            //set its best next info
            predstateinfo->bestnextstate = MDPstate;

            //check the decrease of g-values along the path
            if (predstateinfo->v >= stateinfo->g) {
                SBPL_ERROR("ERROR in ReconstructPath: g-values are non-decreasing\n");
                PrintSearchState(predstateinfo, fDeb);
                throw new SBPL_Exception();
            }

            //transition back
            MDPstate = PredMDPstate;
        }
    }

    return 1;
}

void EPASEPlanner::PrintSearchPath(EPASESearchStateSpace_t* pSearchStateSpace, FILE* fOut)
{
    EPASEState* searchstateinfo;
    CMDPSTATE* state;
    int goalID;
    int PathCost;

    if (bforwardsearch) {
        state = pSearchStateSpace->searchstartstate;
        goalID = pSearchStateSpace->searchgoalstate->StateID;
    }
    else {
        state = pSearchStateSpace->searchgoalstate;
        goalID = pSearchStateSpace->searchstartstate->StateID;
    }
    if (fOut == NULL) fOut = stdout;

    PathCost = ((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;

    SBPL_FPRINTF(fOut, "Printing a path from state %d to the goal state %d\n", state->StateID,
                 pSearchStateSpace->searchgoalstate->StateID);
    SBPL_FPRINTF(fOut, "Path cost = %d:\n", PathCost);

    environment_->PrintState(state->StateID, false, fOut);

    int costFromStart = 0;
    while (state->StateID != goalID) {
        SBPL_FPRINTF(fOut, "state %d ", state->StateID);

        if (state->PlannerSpecificData == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (EPASEState*)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }
        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        int costToGoal = PathCost - costFromStart;
        int transcost = searchstateinfo->g - ((EPASEState*)(searchstateinfo->bestnextstate->PlannerSpecificData))->v;
        if (bforwardsearch) transcost = -transcost;

        costFromStart += transcost;

        SBPL_FPRINTF(fOut, "g=%d-->state %d, h = %d ctg = %d  ", searchstateinfo->g,
                     searchstateinfo->bestnextstate->StateID, searchstateinfo->h, costToGoal);

        state = searchstateinfo->bestnextstate;

        environment_->PrintState(state->StateID, false, fOut);
    }
}

void EPASEPlanner::PrintSearchState(EPASEState* state, FILE* fOut)
{
#if DEBUG
    SBPL_FPRINTF(fOut, "state %d: h=%d g=%u v=%u iterc=%d callnuma=%d expands=%d heapind=%d inconslist=%d\n",
                 state->MDPstate->StateID, state->h, state->g, state->v,
                 state->iterationclosed, state->callnumberaccessed, state->numofexpands, state->heapindex,
                 state->listelem[EPASE_INCONS_LIST_ID] ? 1 : 0);
#else
    SBPL_FPRINTF(fOut, "state %d: h=%d g=%u v=%u iterc=%d callnuma=%d heapind=%d inconslist=%d\n",
                 state->MDPstate->StateID, state->h, state->g, state->v, state->iterationclosed,
                 state->callnumberaccessed, state->heapindex, state->listelem[EPASE_INCONS_LIST_ID] ? 1 : 0);
#endif
    environment_->PrintState(state->MDPstate->StateID, true, fOut);
}

int EPASEPlanner::getHeurValue(EPASESearchStateSpace_t* pSearchStateSpace, int StateID)
{
    CMDPSTATE* MDPstate = GetState(StateID, pSearchStateSpace);
    EPASEState* searchstateinfo = (EPASEState*)MDPstate->PlannerSpecificData;
    return searchstateinfo->h;
}

vector<int> EPASEPlanner::GetSearchPath(EPASESearchStateSpace_t* pSearchStateSpace, int& solcost)
{
    vector<int> SuccIDV;
    vector<int> CostV;
    vector<int> wholePathIds;
    EPASEState* searchstateinfo;
    CMDPSTATE* state = NULL;
    CMDPSTATE* goalstate = NULL;
    CMDPSTATE* startstate = NULL;

    if (bforwardsearch) {
        startstate = pSearchStateSpace->searchstartstate;
        goalstate = pSearchStateSpace->searchgoalstate;

        //reconstruct the path by setting bestnextstate pointers appropriately
        ReconstructPath(pSearchStateSpace);
    }
    else {
        startstate = pSearchStateSpace->searchgoalstate;
        goalstate = pSearchStateSpace->searchstartstate;
    }

    state = startstate;

    wholePathIds.push_back(state->StateID);
    solcost = 0;

    FILE* fOut = stdout;
    if (fOut == NULL) {
        SBPL_ERROR("ERROR: could not open file\n");
        throw new SBPL_Exception();
    }
    while (state->StateID != goalstate->StateID) {
        if (state->PlannerSpecificData == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (EPASEState*)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }
        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        environment_->GetSuccs(state->StateID, &SuccIDV, &CostV);
        int actioncost = INFINITECOST;
        for (int i = 0; i < (int)SuccIDV.size(); i++) {
            if (SuccIDV.at(i) == searchstateinfo->bestnextstate->StateID && CostV.at(i) < actioncost) {
                actioncost = CostV.at(i);
            }
        }
        if (actioncost == INFINITECOST) SBPL_PRINTF("WARNING: actioncost = %d\n", actioncost);

        solcost += actioncost;

        //SBPL_FPRINTF(fDeb, "actioncost=%d between states %d and %d\n",
        //        actioncost, state->StateID, searchstateinfo->bestnextstate->StateID);
        //environment_->PrintState(state->StateID, false, fDeb);
        //environment_->PrintState(searchstateinfo->bestnextstate->StateID, false, fDeb);

#if DEBUG
        EPASEState* nextstateinfo = (EPASEState*)(searchstateinfo->bestnextstate->PlannerSpecificData);
        if (actioncost != abs((int)(searchstateinfo->g - nextstateinfo->g)) && 
            pSearchStateSpace->eps_satisfied <= 1.001)
        {
            SBPL_FPRINTF(fDeb, "ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                         actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            SBPL_ERROR("ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                       actioncost,abs((int)(searchstateinfo->g - nextstateinfo->g)));
            PrintSearchState(searchstateinfo, fDeb);
            PrintSearchState(nextstateinfo, fDeb);
        }
#endif

        state = searchstateinfo->bestnextstate;

        wholePathIds.push_back(state->StateID);
    }

    return wholePathIds;
}

bool EPASEPlanner::Search(EPASESearchStateSpace_t* pSearchStateSpace, vector<int>& pathIds, int & PathCost,
                        bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs)
{
    CKey key;
    timeStart = Clock::now();
    searchexpands = 0;
    num_of_expands_initial_solution = -1;
    double old_repair_time = repair_time;
    if (!use_repair_time)
        repair_time = MaxNumofSecs;

#if DEBUG
    SBPL_FPRINTF(fDeb, "new search call (call number=%d)\n", pSearchStateSpace->callnumber);
#endif

    if (pSearchStateSpace->bReevaluatefvals) {
        // costs have changed or a new goal has been set
        environment_->EnsureHeuristicsUpdated(bforwardsearch);
        Reevaluatehvals(pSearchStateSpace);
    }

    if (pSearchStateSpace->bReinitializeSearchStateSpace) {
        //re-initialize state space
        ReInitializeSearchStateSpace(pSearchStateSpace);
    }

    if (bOptimalSolution) {
        pSearchStateSpace->eps = 1;
        MaxNumofSecs = INFINITECOST;
        repair_time = INFINITECOST;
    }
    else if (bFirstSolution) {
        MaxNumofSecs = INFINITECOST;
        repair_time = INFINITECOST;
    }

    //the main loop of EPASE*
    stats.clear();
    int prevexpands = 0;
    Clock::time_point loop_time;
    while (pSearchStateSpace->eps_satisfied > final_epsilon &&
           GetTimeSince(timeStart) < MaxNumofSecs &&
               (pSearchStateSpace->eps_satisfied == INFINITECOST ||
               GetTimeSince(timeStart) < repair_time ))
    {
        loop_time = Clock::now();
        //decrease eps for all subsequent iterations
        if (fabs(pSearchStateSpace->eps_satisfied - pSearchStateSpace->eps) < ERR_EPS && !bFirstSolution) {
            pSearchStateSpace->eps = pSearchStateSpace->eps - dec_eps;
            if (pSearchStateSpace->eps < final_epsilon)
                pSearchStateSpace->eps = final_epsilon;

            //the priorities need to be updated
            pSearchStateSpace->bReevaluatefvals = true;

            //it will be a new search
            pSearchStateSpace->bNewSearchIteration = true;
        }           

        if (pSearchStateSpace->bNewSearchIteration) {
            pSearchStateSpace->searchiteration++;
            pSearchStateSpace->bNewSearchIteration = false;
            BuildNewOPENList(pSearchStateSpace);
        }

        //re-compute f-values if necessary and reorder the heap
        if (pSearchStateSpace->bReevaluatefvals)
            Reevaluatefvals(pSearchStateSpace);

        //improve or compute path
        if (ImprovePath(pSearchStateSpace, MaxNumofSecs) == 1) {
            pSearchStateSpace->eps_satisfied = pSearchStateSpace->eps;
        }

        //print the solution cost and eps bound
        SBPL_PRINTF("eps=%f expands=%d g(searchgoal)=%d time=%.3f\n", pSearchStateSpace->eps_satisfied,
                    searchexpands - prevexpands,
                    ((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g,
                    GetTimeSince(loop_time));

        if (pSearchStateSpace->eps_satisfied == finitial_eps && pSearchStateSpace->eps == finitial_eps) {
            finitial_eps_planning_time = GetTimeSince(loop_time);
            num_of_expands_initial_solution = searchexpands - prevexpands;
        }

        if (stats.empty() || pSearchStateSpace->eps_satisfied != stats.back().eps) {
            PlannerStats tempStat;
            tempStat.eps = pSearchStateSpace->eps_satisfied;
            tempStat.expands = searchexpands - prevexpands;
            tempStat.time = GetTimeSince(loop_time);
            tempStat.cost = ((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
            stats.push_back(tempStat);
        }

#if DEBUG
        SBPL_FPRINTF(fDeb, "eps=%f expands=%d g(searchgoal)=%d time=%.3f\n", pSearchStateSpace->eps_satisfied,
                     searchexpands - prevexpands,
                     ((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g,
                     GetTimeSince(loop_time));
        PrintSearchState((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData, fDeb);
#endif
        prevexpands = searchexpands;

        //if just the first solution then we are done
        if (bFirstSolution)
            break;

        //no solution exists
        if (((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g == INFINITECOST)
            break;
    }
    repair_time = old_repair_time;

#if DEBUG
    SBPL_FFLUSH(fDeb);
#endif

    PathCost = ((EPASEState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
    MaxMemoryCounter += environment_->StateID2IndexMapping.size() * sizeof(int);

    SBPL_PRINTF("MaxMemoryCounter = %d\n", MaxMemoryCounter);

    int solcost = INFINITECOST;
    bool ret = false;
    if (PathCost == INFINITECOST) {
        SBPL_PRINTF("could not find a solution\n");
        ret = false;
    }
    else {
        SBPL_PRINTF("solution is found\n");
        pathIds = GetSearchPath(pSearchStateSpace, solcost);
        ret = true;
    }

    SBPL_PRINTF("total expands this call = %d, planning time = %.3f secs, solution cost=%d\n",
                searchexpands, GetTimeSince(timeStart), solcost);
    final_eps_planning_time = GetTimeSince(timeStart);
    final_eps = pSearchStateSpace->eps_satisfied;
    //SBPL_FPRINTF(fStat, "%d %d\n", searchexpands, solcost);

    return ret;
}

//-----------------------------Interface function-----------------------------------------------------

int EPASEPlanner::replan(vector<int>* solution_stateIDs_V, ReplanParams params)
{
    int solcost;
    return replan(solution_stateIDs_V, params, &solcost);
}

int EPASEPlanner::replan(vector<int>* solution_stateIDs_V, ReplanParams params, int* solcost)
{
    finitial_eps = params.initial_eps;
    final_epsilon = params.final_eps;
    dec_eps = params.dec_eps;
    bsearchuntilfirstsolution = params.return_first_solution;
    use_repair_time = params.repair_time > 0;
    repair_time = params.repair_time;
    return replan(params.max_time, solution_stateIDs_V, solcost);
}

//returns 1 if found a solution, and 0 otherwise
int EPASEPlanner::replan(double allocated_time_secs, vector<int>* solution_stateIDs_V)
{
    int solcost;

    return replan(allocated_time_secs, solution_stateIDs_V, &solcost);
}

//returns 1 if found a solution, and 0 otherwise
int EPASEPlanner::replan(double allocated_time_secs, vector<int>* solution_stateIDs_V, int* psolcost)
{
    vector<int> pathIds;
    bool bFound = false;
    int PathCost;
    bool bFirstSolution = this->bsearchuntilfirstsolution;
    bool bOptimalSolution = false;
    *psolcost = 0;

    SBPL_PRINTF("planner: replan called (bFirstSol=%d, bOptSol=%d)\n", bFirstSolution, bOptimalSolution);

    //plan
    if (!(bFound = Search(pSearchStateSpace_, pathIds, PathCost,
                          bFirstSolution, bOptimalSolution, allocated_time_secs)))
    {
        SBPL_PRINTF("failed to find a solution\n");
    }

    //copy the solution
    *solution_stateIDs_V = pathIds;
    *psolcost = PathCost;

    return (int)bFound;
}

int EPASEPlanner::set_goal(int goal_stateID)
{
    SBPL_PRINTF("planner: setting goal to %d\n", goal_stateID);
    environment_->PrintState(goal_stateID, true, stdout);

    if (bforwardsearch) {
        if (SetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }
    }
    else {
        if (SetSearchStartState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search start state\n");
            return 0;
        }
    }

    return 1;
}

int EPASEPlanner::set_start(int start_stateID)
{
    SBPL_PRINTF("planner: setting start to %d\n", start_stateID);
    environment_->PrintState(start_stateID, true, stdout);

    if (bforwardsearch) {
        if (SetSearchStartState(start_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search start state\n");
            return 0;
        }
    }
    else {
        if (SetSearchGoalState(start_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }
    }

    return 1;
}

void EPASEPlanner::costs_changed(StateChangeQuery const & stateChange)
{
    pSearchStateSpace_->bReevaluatefvals = true;
    pSearchStateSpace_->bReinitializeSearchStateSpace = true;
}

void EPASEPlanner::costs_changed()
{
    pSearchStateSpace_->bReevaluatefvals = true;
    pSearchStateSpace_->bReinitializeSearchStateSpace = true;
}

int EPASEPlanner::force_planning_from_scratch()
{
    SBPL_PRINTF("planner: forceplanfromscratch set\n");

    pSearchStateSpace_->bReinitializeSearchStateSpace = true;

    return 1;
}

int EPASEPlanner::force_planning_from_scratch_and_free_memory()
{
    SBPL_PRINTF("planner: forceplanfromscratch set\n");
    int start_id = -1;
    int goal_id = -1;
    if (pSearchStateSpace_->searchstartstate)
        start_id = pSearchStateSpace_->searchstartstate->StateID;
    if (pSearchStateSpace_->searchgoalstate)
        goal_id = pSearchStateSpace_->searchgoalstate->StateID;

    if (!bforwardsearch) {
        int temp = start_id;
        start_id = goal_id;
        goal_id = temp;
    }

    DeleteSearchStateSpace(pSearchStateSpace_);
    CreateSearchStateSpace(pSearchStateSpace_);
    InitializeSearchStateSpace(pSearchStateSpace_);
    for (unsigned int i = 0; i < environment_->StateID2IndexMapping.size(); i++)
        for (int j = 0; j < NUMOFINDICES_STATEID2IND; j++)
            environment_->StateID2IndexMapping[i][j] = -1;

    if (start_id >= 0)
        set_start(start_id);
    if (goal_id >= 0)
        set_goal(goal_id);
    return 1;
}

int EPASEPlanner::set_search_mode(bool bSearchUntilFirstSolution)
{
    SBPL_PRINTF("planner: search mode set to %d\n", bSearchUntilFirstSolution);

    bsearchuntilfirstsolution = bSearchUntilFirstSolution;

    return 1;
}

void EPASEPlanner::print_searchpath(FILE* fOut)
{
    PrintSearchPath(pSearchStateSpace_, fOut);
}

//---------------------------------------------------------------------------------------------------------

void EPASEPlanner::get_search_stats(vector<PlannerStats>* s)
{
    s->clear();
    s->reserve(stats.size());
    for (unsigned int i = 0; i < stats.size(); i++) {
        s->push_back(stats[i]);
    }
}

/*bool EPASEPlanner::outOfTime(){
  double time_used = GetTimeSince(timeStart);
  //we are out of time if:
         //we used up the max time limit OR
         //we found some solution and used up the minimum time limit
  if(params.return_first_solution)
    return false;
  if(time_used >= params.max_time)
    printf("out of max time\n");
  if(use_repair_time && eps_satisfied != INFINITECOST && time_used >= params.repair_time)
    printf("used all repair time...\n");
  return time_used >= params.max_time || 
         (use_repair_time && eps_satisfied != INFINITECOST && time_used >= params.repair_time);
}*/