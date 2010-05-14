/*Copyright 2009 Alex Graves

This file is part of RNNLIB.

RNNLIB is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RNNLIB is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RNNLIB.  If not, see <http://www.gnu.org/licenses/>.*/

#ifndef _INCLUDED_Rprop_h  
#define _INCLUDED_Rprop_h

#include <algorithm>
#include "Optimiser.hpp"
#include "DataExporter.hpp"
#include "Helpers.hpp"

using namespace std;

namespace rnnlib {

struct Rprop: public DataExporter, public Optimiser
{
	//data
	ostream& out;
	vector<double> deltas;
	vector<double> prevDerivs;
	double etaChange;
	double etaMin;
	double etaPlus;
	double minDelta;
	double maxDelta;
	double initDelta;
	double prevAvgDelta;
	bool online;
	vector<double>& derivs;
	vector<double>& wts;
	vector<double>& plasts;
	WeightContainer* weightContainer;

	//functions
	Rprop(ostream& o, WeightContainer* wc, DataExportHandler* dEH, bool on = false, const string& name = "optimiser"):
		DataExporter(name, dEH),
		out(o),
		etaChange(0.01),
		etaMin(0.5),
		etaPlus(1.2),
		minDelta(1e-9),
		maxDelta(0.2),
		initDelta(0.01),
		prevAvgDelta(0),
		online(on),
		derivs(wc->derivatives),
		wts(wc->weights),
		plasts(wc->plasticities)
	{
		weightContainer = wc;
		if (online)
		{
			SAVE(prevAvgDelta);
			SAVE(etaPlus);
		}
		build();
	}
	void update_weights()
	{
		assert(wts.size() == derivs.size());
		assert(wts.size() == deltas.size());
		assert(wts.size() == prevDerivs.size());
		loop (int i, indices(wts))
		{
			double deriv = derivs[i];
			double delta = deltas[i];
			double derivTimesPrev =  deriv * prevDerivs[i];
			if(derivTimesPrev > 0)
			{
				deltas[i] = bound(delta * etaPlus, minDelta, maxDelta);
				wts[i] -= sign(deriv) * delta;
				prevDerivs[i] = deriv;
			}
			else if(derivTimesPrev < 0)
			{
				deltas[i] = bound(delta * etaMin, minDelta, maxDelta);
				prevDerivs[i] = 0;
			}
			else
			{
				wts[i] -= sign(deriv) * delta;
				prevDerivs[i] = deriv;
			}
		}
		//use eta adaptations for online training (from Mike Schuster's thesis)
		if (online)	
		{
 			double avgDelta = mean(deltas);
			if (avgDelta > prevAvgDelta)
			{
				etaPlus = std::max (1.0, etaPlus - etaChange);
			}
			else
			{
				etaPlus += etaChange;
			}
			prevAvgDelta = avgDelta;
		}
		if (GlobalVariables::instance().isVerbose())
		{
			PRINT(minmax(wts), out);
			PRINT(minmax(derivs), out);
			PRINT(minmax(deltas), out);
			PRINT(minmax(prevDerivs), out);
		}
	}
	//NOTE must be called after any change to weightContainer
	void build()
	{
		if (deltas.size() != wts.size())
		{
			deltas.resize(wts.size());
			prevDerivs.resize(wts.size());
			fill(deltas, initDelta);
			fill(prevDerivs, 0);
			weightContainer->save_by_conns(deltas, "deltas");
			weightContainer->save_by_conns(prevDerivs, "prevDerivs");
		}
	}
	void print(ostream& out = cout) const
	{
		out << "RPROP" << endl;
		PRINT(online, out);
		if (online)
		{
			PRINT(prevAvgDelta, out);
			PRINT(etaChange, out);
		} 
		PRINT(etaMin, out);
		PRINT(etaPlus, out);
		PRINT(minDelta, out);
		PRINT(maxDelta, out);
		PRINT(initDelta, out);
	}
};

};

#endif
