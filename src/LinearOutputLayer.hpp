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

#ifndef _INCLUDED_LinearOutputLayer_h  
#define _INCLUDED_LinearOutputLayer_h  

#include "OutputLayer.hpp"
#include "Helpers.hpp"

namespace rnnlib {

struct LinearOutputLayer: public OutputLayer
{	
	//functions
	LinearOutputLayer(const string& name, size_t numSeqDims, size_t size, WeightContainer* wc):
		OutputLayer(name, numSeqDims, size, wc)
	{
		display(this->outputErrors, "errors");
		display(this->outputActivations, "activations");
	}
	void feed_forward(const vector<int>& coords)
	{
		copy(this->inputActivations[coords], this->outputActivations[coords]);
	}
	void feed_back(const vector<int>& coords)
	{
		copy(this->outputErrors[coords], this->inputErrors[coords]);
	}	
};

};

#endif
