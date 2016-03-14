#ifndef __ZEROORDER_CONTROLLER_H__
#define __ZEROORDER_CONTROLLER_H__

#include <armadillo>
#include <cassert>

#include "modelstructure.h"
#include "model.h"

namespace fctrl {

	class ZeroOrder_Controller {
	public:
		ModelStructure * ms;

		double delta;
		double eta;
		double rho;

		arma::mat L;
		arma::mat Linv;

		ZeroOrder_Controller(ModelStructure * ms = 0)
			: ms(ms), delta(0.0001), eta(-1), rho(-1)
		{ }

		bool isDefined() const
		{ return (not Linv.is_empty()); }

		void reset()
		{
			L.reset();
			Linv.reset();
		}

		bool is_delta_balanced() const {
			assert(eta >= 0 and rho >= 0);
			if(eta >= 1/10.0 * delta * rho and eta <= 10 * delta * rho)
				return true;

			return false;
		}

		void reconstruct(const arma::mat& L, const arma::mat& Linv) {
			assert(ms != 0);
			assert(ms->x_order == 0 and ms->y_order == 0);

			this->L = L;
			assert(L.n_rows == ms->vActive_OB_IDs.size());
			assert(L.n_cols == ms->vActive_CP_IDs.size());

			this->Linv = Linv;
			assert(Linv.n_rows == ms->vActive_CP_IDs.size());
			assert(Linv.n_cols == ms->vActive_CP_IDs.size() + ms->vActive_OB_IDs.size());
		}

		arma::colvec get_control_input(const LDS_State& state)
		{
			//Assumption: state representation uses ms->get_OB_center() as modelest_offset_y
			assert(ms->flag_subtract_yoperating_in_estimation_for_DC_offset_correction == true
				or ms->flag_subtract_yobjective_in_estimation_for_DC_offset_correction == true);
			arma::colvec ytm1_sc_error = state.sy;
			arma::colvec utm1          = state.sx;

			std::cout << "get_control_input(): ytm1_sc_error=" << ytm1_sc_error << std::endl;
			std::cout << "utm1=" << utm1 << std::endl;

			arma::colvec y_ext = arma::zeros<arma::colvec>(ytm1_sc_error.n_rows + utm1.n_rows);
			y_ext.subvec(0, ytm1_sc_error.n_rows-1) = ytm1_sc_error;

			arma::colvec incr_ut = Linv * y_ext;
			
			eta = squared( arma::norm(L * incr_ut - ytm1_sc_error, 2) );
			rho = squared( arma::norm(incr_ut, 2) );

			arma::colvec ut = incr_ut + utm1;

			for(int j=0; j<(int)ut.n_cols; j++) {
				if(ut(j) < -1)
					ut(j) = -1;
				if(ut(j) > 1)
					ut(j) = 1;
			}

			return ut;
		}
	};
}

#endif //__ZEROORDER_CONTROLLER_H__
