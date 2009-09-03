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

#ifndef _INCLUDED_SteepestDescent_h  
#define _INCLUDED_SteepestDescent_h

#include <algorithm>
#include "Optimiser.hpp"
#include "DataExporter.hpp"

using namespace std;
//extern bool verbose;

namespace rnnlib {

struct SteepestDescent: public DataExporter, public Optimiser
{
	//data
	ostream& out;
	vector<double> deltas;
	double learnRate;
	double momentum;
	vector<double>& wts;
	vector<double>& derivs;
	vector<double>& plasts;

	//functions
	SteepestDescent(ostream& o, double lr = 1e-4, double mom = 0.9, const string& name = "optimiser"):
		DataExporter(name),
		out(o),
		learnRate(lr),
		momentum(mom),
		wts(WeightContainer::instance().weights),
		derivs(WeightContainer::instance().derivatives),
		plasts(WeightContainer::instance().plasticities)
	{
		build();
	}
	void update_weights()
	{
		assert(wts.size() == derivs.size());
		assert(wts.size() == deltas.size());
		loop (int i, indices(wts))
		{
			double delta = plasts[i] * ((momentum * deltas[i]) - (learnRate * derivs[i]));
			deltas[i] = delta;
			wts[i] += delta;
		}
		if (GlobalVariables::instance().isVerbose())
		{
			out << "weight updates:" << endl;
			PRINT(minmax(wts), out);
			PRINT(minmax(derivs), out);
			PRINT(minmax(deltas), out);
		}
	}
	//NOTE must be called after any change to weightContainer
	void build()
	{
		if (deltas.size() != wts.size())
		{		
			deltas.resize(wts.size());
			fill(deltas, 0);
			WeightContainer::instance().save_by_conns(deltas, "deltas");
		}
	}
	void print(ostream& out = cout) const
	{
		out << "steepest descent" << endl;
		PRINT(learnRate, out);
		PRINT(momentum, out);
	}
};

};

#endif
