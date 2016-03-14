#ifndef __STABILITY_COVERAGE_METRICS_H__
#define __STABILITY_COVERAGE_METRICS_H__

#include <map>
#include <ostream>
#include <set>
#include <vector>

#include "debuglog.h"
#include "frameobservation.h"
#include "model.h"

namespace fctrl {

class HistogramBuilder {
public:

    std::map<size_t, double> range_min;
    std::map<size_t, double> range_max;

    HistogramBuilder()
        : range_min(), range_max()
    {}

    void reset()
    {
        range_min.clear();
        range_max.clear();
    }

    // At timestep t, T_e entries corresponds to timesteps t-1 downto t-|T_e|
    std::vector<size_t> construct_histogram(const std::deque<double>& T_e,
                                            int earliest_time_index,
                                            size_t L,
                                            bool resetBins,
                                            double e_min,
                                            double e_max)
    {
#ifdef FC_STABILITYLENGTH_DEBUG
        std::cout << "--- construct_histogram(): |T_e|=" << T_e.size() << " earliest_time_index=" << earliest_time_index
                  << " L=" << L << " resetBins=" << resetBins << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
        if(range_min.count(L) == 0 || resetBins == true) {
            range_min[L] = e_min;
            range_max[L] = e_max;
        }

        const int numBins = 10;

        std::vector<size_t> bins(numBins, 0);

#ifdef FC_STABILITYLENGTH_DEBUG
        if(earliest_time_index >= int(T_e.size()))
             std::cout << "ERROR: earliest_time_index=" << earliest_time_index << " |T_e|=" << T_e.size() << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
        assert(earliest_time_index < int(T_e.size()));
        assert(earliest_time_index >= int(L) - 1);
        double r_min =  range_min[L];
        double r_max =  range_max[L];
        double r_diff = r_max - r_min;
        assert(r_diff >= 0.0);
        if(r_diff > 0) {
          double one_over_r_diff = 1.0/r_diff;
          for(int i=earliest_time_index; i>=earliest_time_index - int(L) + 1; i--) {
              double e_tpp = T_e[i];
              int bin = static_cast<int>((e_tpp - r_min) * one_over_r_diff * numBins);
              if(bin < 0)
                  bin = 0;
              if(bin >= numBins)
                  bin = numBins - 1;
              bins[bin]++;
          }
        }
        else { // r_diff == 0.0
          bins[0]++;
        }

        return bins;
    }

    //Corresponding bins in bins1 and bins2 are assumed to hold frequency counts
    //over the same range of values. bins1 and bins are histograms over sample segments
    //of length L
    static double kolmogorov_smirnov_D(const std::vector<size_t>& bins1,
                                       const std::vector<size_t>& bins2,
                                       size_t L)
    {
        assert(bins1.size() == bins2.size());
        int maxdiff = 0, cdf1 = 0, cdf2 = 0;

        for(size_t bin=0; bin<bins1.size(); bin++) {
            cdf1 += bins1[bin];
            cdf2 += bins2[bin];
            maxdiff = std::max<int>(maxdiff, std::abs(cdf1 - cdf2));
        }
        assert(cdf1 == int(L));
        assert(cdf2 == int(L));
        double Dmetric = static_cast<double>(maxdiff) / L;
#ifdef FC_STABILITYLENGTH_DEBUG
        std::cout << "KS: Dmetric=" << Dmetric << " for L=" << L << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
        return Dmetric;
    }
}; //class HistogramBuilder

class HistogramStatistics {
public:
    enum stability_type {
        unknown, stable, unstable, highlyunstable
    };

    HistogramBuilder hist_builder;

    std::map<size_t, std::vector<size_t> > last_histogram;
    std::map<size_t, size_t>               last_timestamp;

    std::map<size_t, std::deque<double> >  D_seq;
    std::map<size_t, std::deque<size_t> >  tD_seq;
    std::map<size_t, double>               avgD;

    size_t                                 t_prev_bcp;

    HistogramStatistics()
        : hist_builder(), last_histogram(), last_timestamp(), D_seq(), tD_seq(), avgD(), t_prev_bcp(0)
    {}

    void reset()
    {
        hist_builder.reset();

        last_histogram.clear();
        last_timestamp.clear();

        D_seq.clear();
        tD_seq.clear();
        avgD.clear();
        t_prev_bcp = 0;
    }

    stability_type stability(double Dmetric) const {
        assert(0.0 <= Dmetric && Dmetric <= 1.0);
        if(Dmetric <= 0.10)
            return stable;
        if(Dmetric <= 0.50)
            return unstable;
        return highlyunstable;
    }

    stability_type stability_on_avgD(size_t L) const {
        if(avgD.count(L) == 0)
            return unknown;

        std::map<size_t, double>::const_iterator cit = avgD.find(L), ecit = avgD.end();
        assert(cit != ecit);
        double Dmetric = cit->second;
        return stability( Dmetric ); //i.e. stability( avgD[L] )
    }

    // T_e[0] corresponds to timestep t-1, with progressively older samples in T_e[1], T_e[2], ...
    void update_averaged_kolmogorov_smirnov_distances(bool reconstruct,               //inputs
                                                      const std::deque<double>& T_e,
                                                      const std::set<size_t>& Ls_set,
                                                      double e_min,
                                                      double e_max,
                                                      size_t t_oldest_history_sample,
                                                      size_t t, //current timestep
                                                      size_t& return_Ls, //outputs
                                                      size_t& return_t_bcp)
    {
#ifdef FC_STABILITYLENGTH_DEBUG
        std::cout << "update_averaged_kolmogorov_smirnov_distances(): reconstruct=" << reconstruct << " Ls_set=" << Ls_set << " t=" << t << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
        if(reconstruct == true) {
            last_histogram.clear();
            last_timestamp.clear();
            D_seq.clear();
            tD_seq.clear();
            avgD.clear();
        }

        assert(D_seq.size() == tD_seq.size());
        for(std::set<size_t>::const_iterator it = Ls_set.begin(), eit = Ls_set.end();
            it != eit;
            it++)
        {
            size_t L = *it;
            if(T_e.size() >= L) {
                int hist_timeStep; //timestep of the last (newest) sample in the next segment of length L
                if(last_timestamp.count(L) == 0 ||          //either no previous histogram data
                   last_timestamp[L] < t - T_e.size() - 1)  //or contiguity of histograms over history is now broken
                {
                    //use the timestep of the oldest sample in T_e, + (L - 1) to get to the last sample of the next histogram
                    hist_timeStep =  (int(t) - int(T_e.size())) + (int(L) - 1);
                }
                else {
                    hist_timeStep = int(last_timestamp[L]) + L;
                }

#ifdef FC_STABILITYLENGTH_DEBUG
                std::cout << "last_timestamp[" << L << "]=";
                if(last_timestamp.count(L) == 0)
                    std::cout << "undef";
                else
                    std::cout << last_timestamp[L];
                std::cout << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG

                assert(hist_timeStep >= 0);
                size_t numAppends = 0;
                assert(D_seq[L].size() == tD_seq[L].size());
                for(int tp = hist_timeStep; tp <= int(t)-1; tp += L) {
                    bool resetBins = (D_seq.size() == 0);
                    int earliest_time_index = (int(t) - 1 - tp) + int(L) - 1;
                    std::vector<size_t> newhist = hist_builder.construct_histogram(T_e,
                                                                                   earliest_time_index,
                                                                                   L,
                                                                                   resetBins,
                                                                                   e_min,
                                                                                   e_max);
                    if(last_timestamp.count(L) > 0) { //is the preceding histogram available?
                        double newD = hist_builder.kolmogorov_smirnov_D(newhist, last_histogram[L], L);
                        D_seq[L].push_front(newD);
                        tD_seq[L].push_front( last_timestamp[L] );
                        numAppends++;
                    }
                    last_histogram[L] = newhist;
                    last_timestamp[L] = tp;
                }
                if(numAppends > 0) {
                    const size_t min_number_of_samples = 20;
                    assert(D_seq[L].size() == tD_seq[L].size());
                    if(D_seq[L].size() > min_number_of_samples) {
                        size_t drop_index = 0; //0 ==> nothing to drop
                        for(size_t i=min_number_of_samples; i<D_seq.size(); i++) {
                            if(tD_seq[L][i] < t_oldest_history_sample) {
                                drop_index = i;
                                break;
                            }
                        }
                        if(drop_index > 0) {
                            D_seq[L].resize(drop_index);
                            tD_seq[L].resize(drop_index);
                        }
                    }
                    double w = ( D_seq[L].size() == 1 ? 1.0 : std::pow(0.10, 1.0/(D_seq[L].size() - 1)) );
                    double numer = 0.0, denom = 0.0;
                    double w_pow = 1.0;
                    for(size_t i=0; i<D_seq[L].size(); i++) {
                        numer += w_pow * D_seq[L][i];
                        denom += w_pow;
                        w_pow *= w;
                    }
                    assert(denom > 0.0);
                    double new_avgD_L = numer / denom;
                    avgD[L] = new_avgD_L;
                }
            }
        }

        return_Ls = 0; //0 indicates undefined
        // in ascending order of L in Ls_set
        for(std::set<size_t>::const_iterator it = Ls_set.begin(), eit = Ls_set.end();
            it != eit;
            it++)
        {
            size_t L = *it;
            if(stability_on_avgD(L) == stable) {
                return_Ls = L;
                break;
            }
        }

        return_t_bcp = 0; //0 indicates no behavior change point found
        // in ascending order of L in Ls_set
        if(return_Ls > 0) { //stability established
            assert(D_seq.count(return_Ls) > 0);
                //by definition (return_Ls != 0 ==> avgD.count(return_Ls) == 1 ==> D_seq.count(return_Ls) == 1)
 
            assert(D_seq[return_Ls].size() == tD_seq[return_Ls].size());
            for(size_t i=0; i<D_seq[return_Ls].size(); i++) {
                if(tD_seq[return_Ls][i] < t_oldest_history_sample || tD_seq[return_Ls][i] <= t_prev_bcp)
                    break;
                if(stability( D_seq[return_Ls][i] ) == highlyunstable) {
                    return_t_bcp = tD_seq[return_Ls][i];
                    t_prev_bcp = return_t_bcp;
                    break;
                }
            }
        }
    }
}; //class HistogramStatistics

class StabilityLengthMetrics {
public:
    HistogramStatistics hist_stats;

    Model*                       p_M;
    size_t                       W;
    size_t                       maxorder;

    std::deque<FrameObservation> T_xy;

    LDS_State                    last_state;
    std::deque<double>           T_e;

    std::set<size_t>             Ls_set;
    double                       e_min;
    double                       e_max;
    bool                         reconstruct;

    StabilityLengthMetrics(Model* p_M)
        : hist_stats(), p_M(p_M), W(0), T_xy(), T_e(), Ls_set(), e_min(-1.0), e_max(-1.0), reconstruct(true)
    {}

    void reset()
    {
        hist_stats.reset();

        assert(p_M != 0);

        W = p_M->ms->active_W;
        maxorder = std::max(p_M->ms->x_order, p_M->ms->y_order);

        T_xy.clear();
        last_state = LDS_State();
        T_e.clear();
        Ls_set.clear();
        e_min = -1.0; //< 0 indicates no value held
        e_max = -1.0; //< 0 indicates no value held
        reconstruct = true;

        Ls_set.insert(W);
    }

    size_t round_length_to_W(size_t Lval) const {
        assert(W != 0);
        size_t Lval_mod_W = Lval % W;
        size_t rounded_Lval = Lval_mod_W < W/2 ? Lval - Lval_mod_W : Lval + (W - Lval_mod_W);
        assert(rounded_Lval % W == 0);
        if(rounded_Lval == 0)
            rounded_Lval = W;
        return rounded_Lval;
    }

    void update_stability_length(const FrameObservation& fo,           //inputs
                                 bool newM,
                                 size_t                  t_oldest_history_sample,
                                 size_t&                 return_Ls,    //outputs
                                 size_t&                 return_t_bcp)
    {
        assert(p_M != 0);
        assert(W != 0);
        assert(maxorder != 0);
        assert(Ls_set.size() > 0);

#ifdef FC_STABILITYLENGTH_DEBUG
        std::cout << "update_stability_length(): |T_xy|=" << T_xy.size() << " |T_e|=" << T_e.size()
                  << " erange=[" << e_min << ", " << e_max << "] reconstruct=" << reconstruct
                  << " (x,y,fn)=(" << fo.vCP_observations << ", " << fo.vOB_observations << ", " << fo.frame_number << ") "
                  << " newM=" << newM << " M=" << p_M->isDefined() << " t_oldest_history_sample=" << t_oldest_history_sample
                  << " Ls_set=" << Ls_set << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG

        T_xy.push_front(fo);
        if(newM == true || p_M->isDefined() == false) {
            T_e.clear();
            reconstruct = true;
#ifdef FC_STABILITYLENGTH_DEBUG
            std::cout << "set: T_e.clear(), reconstruct=1" << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
        }

        if(p_M->isDefined() == false) {
            return_Ls = 0;    //indicates undefined
            return_t_bcp = 0; //indicates no behavior change point found

            return;
        }

        if(T_xy.size() >= maxorder){ //i.e. model can be applied to compute at least one output
            arma::colvec CP_N     = p_M->ms->get_CP_N();
            arma::colvec OB_delta = p_M->ms->get_OB_delta();
            arma::colvec OB_importance = p_M->ms->get_OB_importance();

            double prev_e_min = e_min;
            double prev_e_max = e_max;
            if(T_e.size() == 0) {
                //state needs to be setup first
                last_state = LDS_State(p_M->lp, p_M->model_Y_operating);
                int k=int(T_xy.size())-1;
                while(last_state.isStateFullySetup() == false) {
                    assert(k >= 0);
                    arma::colvec x = get_arma_colvec_for_std( T_xy[k].vCP_observations ) / CP_N;
                    arma::colvec y = get_arma_colvec_for_std( T_xy[k].vOB_observations ) / OB_delta;
                    last_state.state_transition(x,y);
                    k--;
                }
                while(k >= 0) {
                    arma::colvec x = get_arma_colvec_for_std( T_xy[k].vCP_observations ) / CP_N;
                    arma::colvec y = get_arma_colvec_for_std( T_xy[k].vOB_observations ) / OB_delta;

                    arma::colvec y_predicted = p_M->lp.get_prediction(last_state, x);
                    arma::colvec importance_error = (y_predicted - y) % OB_importance;
                    double e = squared( arma::norm(importance_error, 2) );

#ifdef FC_STABILITYLENGTH_DEBUG
                    std::cout << "e=" << e << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
                    T_e.push_front(e);
                    
                    if(e_min < 0.0 || e_min > e)
                        e_min = e;
                    if(e_max < e)
                        e_max = e;
                    
                    last_state.state_transition(x,y);
                    k--;
                }

                assert(T_e.size() > 0); //at least one sample should have been added
            }
            else { // T_e.size() > 0 ==> just need to update latest sample (b/c this function is invoked every frame)
                assert(last_state.isStateFullySetup() == true);

                const int k = 0;
                arma::colvec x = get_arma_colvec_for_std( T_xy[k].vCP_observations ) / CP_N;
                arma::colvec y = get_arma_colvec_for_std( T_xy[k].vOB_observations ) / OB_delta;

                arma::colvec y_predicted = p_M->lp.get_prediction(last_state, x);
                arma::colvec importance_error = (y_predicted - y) % OB_importance;
                double e = squared( arma::norm(importance_error, 2) );

#ifdef FC_STABILITYLENGTH_DEBUG
                std::cout << "e=" << e << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
                T_e.push_front(e);
                
                if(e_min < 0.0 || e_min > e)
                    e_min = e;
                if(e_max < e)
                    e_max = e;
                
                last_state.state_transition(x,y);
            }
            assert(e_max >= e_min);
            if( (e_max - e_min) > 1.10 * (prev_e_max - prev_e_min) ) {
                reconstruct = true;
#ifdef FC_STABILITYLENGTH_DEBUG
                std::cout << "set: reconstruct=1" << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG
            }
        }

        assert(T_xy.size() == T_e.size() + maxorder);

        size_t Ls_max = *(Ls_set.rbegin());
        std::set<size_t> removeSet, addSet;
        do {
            assert(Ls_set.size() > 0);

            hist_stats.update_averaged_kolmogorov_smirnov_distances(reconstruct,
                                                                    T_e,
                                                                    Ls_set,
                                                                    e_min,
                                                                    e_max,
                                                                    t_oldest_history_sample,
                                                                    fo.frame_number,
                                                                    return_Ls,     //outputs
                                                                    return_t_bcp);

            reconstruct = false;
            removeSet.clear();
            addSet.clear();
            size_t L1 = 0, L2 = 0; //invalid values
            Ls_max = *(Ls_set.rbegin());

            //must be in ascending order of L3 in Ls_set
            for(std::set<size_t>::const_iterator it = Ls_set.begin(), eit = Ls_set.end();
                it != eit;
                it++)
            {
                size_t L3 = *it;
                assert(hist_stats.avgD.count(0) == 0);
                    //invalid values must not be present in avgD

                HistogramStatistics::stability_type sL1 = hist_stats.stability_on_avgD(L1);
                HistogramStatistics::stability_type sL2 = hist_stats.stability_on_avgD(L2);
                HistogramStatistics::stability_type sL3 = hist_stats.stability_on_avgD(L3);

                if(sL1 == sL2 && sL2 == sL3 && sL3 != HistogramStatistics::unknown) {
                    removeSet.insert(L2);
                }

                if(sL2 != HistogramStatistics::unknown && sL3 != HistogramStatistics::unknown
                        && sL2 != sL3)
                {
                    size_t L = round_length_to_W( (L2 + L3) / 2);
                    if(L != L2 && L != L3) {
                        assert(L != 0);
                        addSet.insert(L);
                    }
                }
                L1 = L2;
                L2 = L3;
            }

            if(return_Ls == 0 && T_e.size() > 2 * Ls_max)
                addSet.insert(2 * Ls_max);

#ifdef FC_STABILITYLENGTH_DEBUG
            std::cout << "return_Ls=" << return_Ls << " |T_e|=" << T_e.size() << " Ls_max=" << Ls_max << " removeSet=" << removeSet << " addSet=" << addSet << std::endl;
#endif //FC_STABILITYLENGTH_DEBUG

            for(std::set<size_t>::const_iterator it = removeSet.begin(), eit = removeSet.end();
                it != eit;
                it++)
            {
                size_t L = *it;
                assert(L >= W);
                Ls_set.erase( L );
            }

            for(std::set<size_t>::const_iterator it = addSet.begin(), eit = addSet.end();
                it != eit;
                it++)
            {
                size_t L = *it;
                assert(L >= W);
                Ls_set.insert( L );
            }


        } while( removeSet.size() > 0 || addSet.size() > 0 );

        if(T_e.size() > 2 * Ls_max) {
            size_t numSamplesToDrop = T_e.size() - 2 * Ls_max;

            T_xy.resize(T_xy.size() - numSamplesToDrop);
            T_e.resize(T_e.size() - numSamplesToDrop);
        }
    }

    HistogramStatistics::stability_type interpolated_stability(size_t L) const {
        size_t Lp = round_length_to_W(L);
        if(Ls_set.count(Lp) > 0 && hist_stats.stability_on_avgD(Lp) == HistogramStatistics::stable) {
            return HistogramStatistics::stable;
        }

        size_t L1 = 0;
        //explore Ls_set in ascending order
        for(std::set<size_t>::const_iterator cit = Ls_set.begin(), ecit = Ls_set.end();
            cit != ecit;
            cit++)
        {
            size_t L2 = *cit;
            if(L1 != 0 && L1 < Lp && Lp < L2) {
                HistogramStatistics::stability_type s1 = hist_stats.stability_on_avgD(L1);
                HistogramStatistics::stability_type s2 = hist_stats.stability_on_avgD(L2);
                if(s1 == s2)
                    return s1;
            }
            if(L2 != Lp) //in case Lp is present in Ls_set, we want to use its pred and succ as L1 and L2
                L1 = L2;
            if(L1 >= Lp)
                break;
        }
        return HistogramStatistics::unknown;
    }

    void add_stability_candidate(size_t L) {
        size_t Lp = round_length_to_W(L);
        if(Ls_set.count(Lp) == 0)
            Ls_set.insert(Lp);
    }

}; // class StabilityLengthMetrics


struct LengthTimestampTuple {
    size_t length;
    size_t timestamp;

    LengthTimestampTuple(size_t length, size_t timestamp)
        : length(length), timestamp(timestamp)
    {}
};

inline
bool operator<(LengthTimestampTuple ltt1, LengthTimestampTuple ltt2) {
    return ltt1.length < ltt2.length || (ltt1.length == ltt2.length && ltt1.timestamp < ltt2.timestamp);
}

class CoverageLengthMetrics {
public:
    Model*                       p_M;
    size_t                       W;
    size_t                       dims_x;

    std::set<LengthTimestampTuple>  Lc_seq;  //ordered by length (ascending)
    std::list<LengthTimestampTuple> Lc_list; //ordered by timestamps (older to later)
    size_t Lc; // = 0 ==> undefined

    CoverageLengthMetrics(Model* p_M)
        : p_M(p_M), W(0), dims_x(0), Lc_seq(), Lc_list(), Lc(0)
    {}

    void reset() {
        assert(p_M != 0);

        W = p_M->ms->active_W;
        dims_x = p_M->ms->get_CP_N().n_rows;

        Lc_seq.clear();
        Lc_list.clear();

        Lc = W;
    }

    size_t update_coverage_length(bool prev_frame_ended_PFE,
                                  double kappa,
                                  size_t hist_length,
                                  size_t t_bcp,
                                  size_t frame_number)
    {
        if(Lc_seq.size() == 0) {
            Lc_seq.insert( LengthTimestampTuple(W, frame_number) );
            Lc_list.push_back( LengthTimestampTuple(W, frame_number) );
        }

#ifdef FC_COVERAGELENGTH_DEBUG
        std::cout << "update_coverage_length(): Lc_seq=" << Lc_seq << " Lc_list=" << Lc_list << " Lc=" << Lc << std::endl;
#endif //FC_COVERAGELENGTH_DEBUG
        assert(Lc_seq.size() == Lc_list.size());


        //Drop samples with timestamps <= t_bcp
        if(t_bcp != 0) {
            std::list<LengthTimestampTuple>::iterator it = Lc_list.begin();
            while(it != Lc_list.end() && it->timestamp <= t_bcp) {
                std::list<LengthTimestampTuple>::iterator del_it = it;
                it++;
                Lc_seq.erase( *del_it );
                Lc_list.erase( del_it );
            }

            if(Lc_seq.size() == 0) {
                Lc_seq.insert( LengthTimestampTuple(W, frame_number) );
                Lc_list.push_back( LengthTimestampTuple(W, frame_number) );
            }
        }

        //Add coverage length sample at end of PFE
        if(prev_frame_ended_PFE == true) {
            size_t L = 0;
            if(kappa < 1)
                L = hist_length + (1 - kappa) * 2 * dims_x;
            else
                L = std::max( (hist_length >= 2 ? hist_length - 2 : 0), W );
            Lc_seq.insert( LengthTimestampTuple(L, frame_number) );
            Lc_list.push_back( LengthTimestampTuple(L, frame_number) );
        }

        //Drop oldest samples until num_to_retain left
        const int num_to_retain = 20;
        int num_left_to_drop = int(Lc_seq.size()) - num_to_retain;
        const int num_to_drop = num_left_to_drop;
        if(num_left_to_drop > 0) {
            std::list<LengthTimestampTuple>::iterator it = Lc_list.begin();
            while(num_left_to_drop-- > 0) {
                assert(it != Lc_list.end());
                std::list<LengthTimestampTuple>::iterator del_it = it;
                it++;
                Lc_seq.erase( *del_it );
                Lc_list.erase( del_it );
            }
        }

        //Update Lc if end of PFE or samples dropped:
        // Lc = average the 50% to 75% median lengths, once enough samples exist
        if(prev_frame_ended_PFE == true || num_to_drop > 0) {
            const size_t num_to_skip_initially = Lc_seq.size() > 4 ? Lc_seq.size() / 2 : Lc_seq.size() / 4;
            const size_t num_to_use_in_middle  = Lc_seq.size() > 8 ? (Lc_seq.size() / 4 + 1) : Lc_seq.size() - num_to_skip_initially;
            assert(num_to_skip_initially + num_to_use_in_middle <= Lc_seq.size());
            assert(num_to_use_in_middle > 0);

            std::set<LengthTimestampTuple>::iterator it = Lc_seq.begin();
            std::advance(it, num_to_skip_initially);

            size_t Lc_sum = 0;
            for(size_t i=0; i<num_to_use_in_middle; i++) {
                Lc_sum += it->length;
                it++;
            }

            Lc = Lc_sum / num_to_use_in_middle;
        }

        return Lc;
    }
};

} //namespace fctrl

inline
std::ostream& operator<<(std::ostream& os, const fctrl::LengthTimestampTuple& ltt)
{
    os << "(" << ltt.length << ", " << ltt.timestamp << ") ";
    return os;
}


#endif //__STABILITY_COVERAGE_METRICS_H__
