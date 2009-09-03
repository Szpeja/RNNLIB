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

#ifndef _INCLUDED_Helpers_h  
#define _INCLUDED_Helpers_h  

#include <boost/date_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/array.hpp>
#include <boost/timer.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/range.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#include <math.h>
#include <numeric>
#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <set>
#include <algorithm>
#include <iterator>
#include <map>
#include <assert.h>
#include "list_of.hpp"

using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace boost::posix_time;
using namespace boost::gregorian;

namespace rnnlib {

#define loop BOOST_FOREACH
#define loop_back BOOST_REVERSE_FOREACH

typedef std::vector<size_t>::const_iterator VSTCI;
typedef std::vector<double>::iterator VDI;
typedef std::vector<double>::const_iterator VDCI;
typedef std::vector<double>::reverse_iterator VDRI;
typedef std::string::iterator SI;
typedef std::string::const_iterator SCI;
typedef std::vector<int>::iterator VII;
typedef std::vector<std::string>::iterator VSI;
typedef std::vector<std::string>::const_iterator VSCI;
typedef std::vector<int>::reverse_iterator VIRI;
typedef std::vector<std::vector<int> >::reverse_iterator VVIRI;
typedef std::vector<int>::const_iterator VICI;
typedef std::vector<bool>::iterator VBI;
typedef std::vector<float>::iterator VFI;
typedef std::vector<float>::const_iterator VFCI;
typedef std::vector<std::vector<double> >::iterator VVDI;
typedef std::vector<std::vector<double> >::const_iterator VVDCI;
typedef std::vector<std::vector<int> >::iterator VVII;
typedef std::vector<std::vector<int> >::const_iterator VVICI;
typedef std::vector<unsigned int>::iterator VUII;
typedef std::vector<std::vector<float> >::iterator VVFI;
typedef std::map<std::string, std::string>::iterator MSSI;
typedef std::map<std::string, std::string>::const_iterator MSSCI;
typedef std::map<std::string, double>::iterator MSDI;
typedef std::map<std::string, double>::const_iterator MSDCI;
typedef std::map<std::string, std::pair<int,double> >::iterator MSPIDI;
typedef std::map<std::string, std::pair<int,double> >::const_iterator MSPIDCI;
typedef std::vector< std::map<std::string, std::pair<int,double> > >::const_iterator VMSDCI;
typedef std::vector<std::map<std::string, std::pair<int,double> > >::iterator VMSDI;
typedef std::vector<std::map<std::string, std::pair<int,double> > >::reverse_iterator VMSDRI;
typedef std::map<std::string, int>::iterator MSII;
typedef std::map<std::string, int>::const_iterator MSICI;
typedef std::map<int, int>::iterator MIII;
typedef std::map<int, int>::const_iterator MIICI;
typedef std::vector<std::vector<int> >::const_reverse_iterator VVIRCI;
typedef std::vector<int>::const_reverse_iterator VIRCI;
typedef std::vector<const float*>::const_iterator VPCFCI;
typedef std::vector<const float*>::iterator VPCFI;
typedef std::vector<const float*>::const_reverse_iterator VPCFCRI;
typedef std::vector<bool>::const_iterator VBCI;
typedef std::vector<bool>::iterator VBI;
typedef std::map <std::string, std::pair<double, int> >::iterator MCSPDII;
typedef std::map <std::string, std::pair<double, int> >::const_iterator MCSPDICI;
typedef bimap<int, std::string>::left_const_iterator BMISLCI;
typedef bimap<int, std::string>::right_const_iterator BMISRCI;
typedef bimap<int, std::string>::relation BMISR;
typedef std::pair<std::string, double> PSD;
typedef std::pair<int, int> PII;
typedef std::pair<const std::string, double> PCSD;
typedef std::pair<std::string, int> PSI;
typedef std::pair<std::string, std::string> PSS;
typedef const tuple<double&, double&, double&, double&>& TDDDD;
typedef const tuple<double&, double&, double&, double&, double&>& TDDDDD;
typedef const tuple<double&, double&, double&>& TDDD;
typedef const tuple<double&, double&, int&>& TDDI;
typedef const tuple<double&, double&, float&>& TDDF;
typedef const tuple<double&, double&, float>& TDDCF;
typedef const tuple<double&, double&>& TDD;
typedef const tuple<std::string, int>& TSI;
typedef const tuple<int, int>& TII;
typedef const tuple<int, set<int>&>& TISETI;

//global variables
// static const double doubleMax = numeric_limits<double>::max();
// static const double doubleMin = numeric_limits<double>::min();
// static const double infinity = numeric_limits<double>::infinity();
// static bool runningGradTest = false;
// static bool verbose = false;
//static ostream& COUT = cout;

//singleton added by Sergio - only one instance of this class is possible
class GlobalVariables { 
 private:
	double doubleMax;
	double doubleMin;
	double infinity;
	bool runningGradTest;
	bool verbose;
	double expLimit;
	double negExpLimit;
	double logZero;
 public:
	static GlobalVariables& instance() 
	{
		static GlobalVariables inst;
		return inst;
	}
	double getDoubleMax () { return doubleMax; }
	double getDoubleMin () { return doubleMin; }
	double getInfinity () { return infinity; }
	double getExpLimit () { return expLimit; }
	double getNegExpLimit () { return negExpLimit; }
	double getLogZero () { return logZero; }
	bool isRunningGradTest () { return runningGradTest; }
	bool isVerbose () { return verbose; }
	void setRunningGradTest (bool val) { runningGradTest = val; }
	void setVerbose (bool val) { verbose = val; }
	
 protected:
	GlobalVariables() {
		doubleMax = numeric_limits<double>::max();
		doubleMin = numeric_limits<double>::min();
		infinity = numeric_limits<double>::infinity();
		expLimit = log(numeric_limits<double>::max());
		negExpLimit = log(numeric_limits<double>::min());
		logZero = -infinity;
		verbose = false;
		runningGradTest = false;
	}
	virtual ~GlobalVariables() {}
	GlobalVariables(const GlobalVariables&);                 // Prevent copy-construction
	GlobalVariables& operator=(const GlobalVariables&);      // Prevent assignment

};

#define PRINT(x, o) ((o) << boolalpha << #x " = " << (x) << endl)
#define PRINTN(x, o) (o) << boolalpha << #x ":" << endl; print_range((o), (x), std::string("\n")); (o) << endl
#define PRT(x) PRINT(x, cout)
#define PRTN(x) PRINTN(x, cout)
#define PRINTR(x, o) (o) << boolalpha << #x " = "; print_range((o), (x)); (o) << endl
#define PRTR(x) PRINTR(x, cout)
#define check(condition, str)  if(!(condition)) {cout << "ERRROR: " << (str) << endl; assert((condition));}

//MISC FUNCTIONS
static bool warn (bool condition, ostream& out, const std::string& str)
{
	if (!condition)
	{
		out << "WARNING: " << str << endl;
	}
	return condition;
}
static void print_time(double totalSeconds, ostream& out = cout, bool abbrv = false)
{
	int wholeSeconds = floor(totalSeconds);
	int seconds = wholeSeconds % 60;
	int totalMinutes = wholeSeconds / 60;
	int minutes = totalMinutes % 60;
	int totalHours = totalMinutes / 60;
	int hours = totalHours % 24;
	int totalDays = totalHours / 24;
	int days = totalDays % 365;
	if (days)
	{
		out << days << " day";
		if (days > 1)
		{
			out << "s";
		}
		out << " ";
	}
	if (hours)
	{
		out << hours << (abbrv ? " hr" : " hour");
		if (hours > 1)
		{
			out << "s";
		}
		out << " ";
	}
	if (minutes)
	{
		out << minutes << (abbrv ? " min" : " minute");
		if (minutes > 1)
		{
			out << "s";
		}
		out << " ";
	}
	out << totalSeconds - wholeSeconds + seconds << (abbrv ? " secs" : " seconds");
}
static std::string time_stamp(const std::string& format = "%Y.%m.%d-%H.%M.%S%F%Q")
{
	time_facet* timef = new time_facet(format.c_str());
	std::stringstream ss;
	ss.imbue(locale(ss.getloc(), timef));
	ss << microsec_clock::local_time();
	return ss.str();
}
static void mark()
{
	static int num = 0;
	cout << "MARK " << num << endl;
	++num;
}
template<class T> static T squared(const T& t)
{
	return t*t;
}
template<class T> static int sign(const T& t)
{
	if (t < 0)
	{
		return -1;
	}
	else if (t > 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
template <class T> static T bound (const T& v, const T& minVal, const T& maxVal)
{
	return min(max(minVal, v), maxVal);
}
//CAST OPERATIONS
template<class T> static std::string str(const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
//	return lexical_cast<std::string>(t);
}
template<class T> static double dbl(const T& t)
{
	return lexical_cast<double>(t);
}
template<class T> static double flt(const T& t)
{
	return lexical_cast<float>(t);
}
template<class T> static double integer(const T& t)
{
	return lexical_cast<int>(t);
}
//GENERIC RANGE OPERATIONS
template <class R> static size_t count_adjacent(const R& r)
{
	size_t count = 0;
	for (typename range_iterator<R>::type it = boost::begin(r); it != boost::end(r); it = adjacent_find(it, boost::end(r)))
	{
		++count;	
	}
	return count;
}
template <class R1, class R2> static std::pair<zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type> >,
											zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type> > > 
zip(R1& r1, R2& r2)
{
	size_t size = range_min_size(r1, r2);
	return make_pair(make_zip_iterator(make_tuple(boost::begin(r1), boost::begin(r2))), 
					 make_zip_iterator(make_tuple(boost::end(r1) - (boost::size(r1) - size), boost::end(r2) - (boost::size(r2) - size))));
}
template <class R1, class R2, class R3> static std::pair<zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type, typename range_iterator<R3>::type> >,
													zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type, typename range_iterator<R3>::type> > > 
zip(R1& r1, R2& r2, R3& r3)
{
	size_t size = range_min_size(r1, r2, r3);
	return make_pair(make_zip_iterator(make_tuple(boost::begin(r1), boost::begin(r2), boost::begin(r3))), 
					 make_zip_iterator(make_tuple(boost::end(r1) - (boost::size(r1) - size), boost::end(r2) - (boost::size(r2) - size), boost::end(r3) - (boost::size(r3) - size))));
}
template <class R1, class R2, class R3, class R4> static std::pair<zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type, typename range_iterator<R3>::type, typename range_iterator<R4>::type> >,
zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type, typename range_iterator<R3>::type, typename range_iterator<R4>::type> > > 
zip(R1& r1, R2& r2, R3& r3, R4& r4)
{
	size_t size = range_min_size(r1, r2, r3, r4);
	return make_pair(make_zip_iterator(make_tuple(boost::begin(r1), boost::begin(r2), boost::begin(r3), boost::begin(r4))), 
					 make_zip_iterator(make_tuple(boost::end(r1) - (boost::size(r1) - size), boost::end(r2) - (boost::size(r2) - size), boost::end(r3) - (boost::size(r3) - size), boost::end(r4) - (boost::size(r4) - size))));
}
template <class R1, class R2, class R3, class R4, class R5> static std::pair<zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type, typename range_iterator<R3>::type, typename range_iterator<R4>::type, typename range_iterator<R5>::type> >,
zip_iterator<tuple<typename range_iterator<R1>::type, typename range_iterator<R2>::type, typename range_iterator<R3>::type, typename range_iterator<R4>::type, typename range_iterator<R5>::type> > > 
zip(R1& r1, R2& r2, R3& r3, R4& r4, R5& r5)
{
	size_t size = range_min_size(r1, r2, r3, r4, r5);
	return make_pair(make_zip_iterator(make_tuple(boost::begin(r1), boost::begin(r2), boost::begin(r3), boost::begin(r4), boost::begin(r5))), 
					 make_zip_iterator(make_tuple(boost::end(r1) - (boost::size(r1) - size), boost::end(r2) - (boost::size(r2) - size), boost::end(r3) - (boost::size(r3) - size), boost::end(r4) - (boost::size(r4) - size), boost::end(r5) - (boost::size(r5) - size))));
}
template <class R> static std::pair<counting_iterator<typename range_difference<R>::type>, counting_iterator<typename range_difference<R>::type> > indices(const R& r)
{
	return range(boost::size(r));
}
template <class R> static std::pair<zip_iterator<tuple<typename range_iterator<R>::type, counting_iterator<typename range_difference<R>::type> > >,
								zip_iterator<tuple<typename range_iterator<R>::type, counting_iterator<typename range_difference<R>::type> > > >
enumerate(R& r)
{
	return make_pair(make_zip_iterator(make_tuple(boost::begin(r), counting_iterator<typename range_difference<R>::type>(0))), 
					make_zip_iterator(make_tuple(boost::end(r), counting_iterator<typename range_difference<R>::type>(boost::size(r)))));
}
template <class T> static std::pair<counting_iterator<T>, counting_iterator<T> > range(const T& t)
{
	return make_pair(counting_iterator<T>(0), counting_iterator<T>(t));
}
template <class T> static std::pair<counting_iterator<T>, counting_iterator<T> > range(const T& t1, const T& t2)
{
	return make_pair(counting_iterator<T>(t1), counting_iterator<T>(t2));
}
template <class R1, class R2, class F> static typename range_iterator<R2>::type transform(const R1& r1, R2& r2, F f)
{
	return std::transform(boost::begin(r1), boost::end(r1), boost::begin(r2), f);
}
template <class R> static typename range_value<R>::type& nth_last(R& r, size_t n = 1)
{
	check(n > 0 && n <= boost::size(r), "nth_last called with n = " + str(n) + " for range of size " + (str(boost::size(r))));
	return *(boost::end(r) - n); 
}
template <class R> size_t last_index(R& r)
{
	return (boost::size(r) - 1); 
}
template <class R, class UnaryFunction> static UnaryFunction for_each(R& r, UnaryFunction f)
{
	return for_each(boost::begin(r), boost::end(r), f); 
}
template <class R, class T> static bool in(const R& r, const T& t)
{
	return find(boost::begin(r), boost::end(r), t) != boost::end(r);
}
template <class R, class T> static size_t index(const R& r, const T& t)
{
	return distance(boost::begin(r), find(boost::begin(r), boost::end(r), t));
}
template <class R> static void reverse(R& r)
{
	reverse(boost::begin(r), boost::end(r));
}
template <class R> static void sort(R& r)
{
	sort(boost::begin(r), boost::end(r));
}
template <class R> static void reverse_sort(R& r)
{
	sort(boost::rbegin(r), boost::rend(r));
}
template <class R> std::pair<typename range_value<R>::type, typename range_value<R>::type> minmax(const R& r)
{
	std::pair<typename range_const_iterator<R>::type, typename range_const_iterator<R>::type> p = minmax_element(boost::begin(r), boost::end(r));
	return make_pair(*p.first, *p.second); 
}
template <class R> static void bound_range (R& r, const typename boost::range_value<R>::type& minVal, const typename boost::range_value<R>::type& maxVal)
{
	for (typename range_iterator<R>::type it = boost::begin(r); it != boost::end(r); ++it) 
	{
		*it = bound(*it, minVal, maxVal);
	}
}
template <class R1, class R2> typename boost::range_value<R1>::type euclidean_squared(const R1& r1, const R2& r2)
{
	typename range_const_iterator<R2>::type b2 = boost::begin(r2); 
	typename range_const_iterator<R1>::type e = boost::end(r1);
	typename boost::range_value<R1>::type d = 0;
	for (typename range_const_iterator<R1>::type b1 = boost::begin(r1); b1 != e; ++b1, ++b2)
	{
		typename boost::range_value<R1>::type diff = *b1-*b2;
		d += diff * diff;
	}
	return d;
}
template<class R> static void range_negate
		(R& r)
{
	transform(boost::begin(r), boost::end(r), boost::begin(r), negate<typename boost::range_value<R>::type>());
}
template<class R> static void fill (R& r, const typename boost::range_value<R>::type& v)
{
	fill(boost::begin(r), boost::end(r), v);
}
template<class R> static size_t count(const R& r, const typename boost::range_value<R>::type& v)
{
	return count(boost::begin(r), boost::end(r), v);
}
template<class R1, class R2> static void copy(const R1& source, R2& dest)
{
	assert(boost::size(dest) >= boost::size(source));
	std::copy(boost::begin(source), boost::end(source), boost::begin(dest));
}
template<class R1, class R2> static void reverse_copy(const R1& source, R2& dest)
{
	reverse_copy(boost::begin(source), boost::end(source), boost::begin(dest));
}
template <class R> static std::vector<typename boost::range_value<R>::type>& flip(const R& r)
{
	static std::vector<typename boost::range_value<R>::type> v;
	v.resize(boost::size(r));
	reverse_copy(r, v);
	return v;
}
template<class R1, class R2> static bool equal(const R1& source, R2& dest)
{
	return ((boost::size(source) == boost::size(dest)) && std::equal(boost::begin(source), boost::end(source), boost::begin(dest)));
}
template<class R> static void shuffle (R& r)
{
	random_shuffle(boost::begin(r), boost::end(r));
}
template <class R> static typename range_value<R>::type max(const R& r)
{
	return *max_element(boost::begin(r), boost::end(r));
}
template <class C, class Tr, class R> static void print_range(basic_ostream<C, Tr>& out, const R& r, const basic_string<C, Tr>& delim = " ")
{
	typename range_const_iterator<R>::type b = boost::begin(r); 
	typename range_const_iterator<R>::type e = boost::end(r);
	if (b != e) 
	{ 
		out << *b;
		while (++b != e) 
		{
			out << delim << *b; 
		}
	}
}
template <class C, class Tr, class R> static basic_ostream<C, Tr>& operator <<(basic_ostream<C, Tr>& out, const R& r)
{
	print_range(out, r);
	return out;
}
template <class C, class Tr, class R> static basic_istream<C, Tr>& operator >>(basic_istream<C, Tr>& in, R& r)
{
	typename range_iterator<R>::type b = boost::begin(r); 
	typename range_iterator<R>::type e = boost::end(r);
	for (; b != e; ++b)
	{
		in >> *b; 
	}
	return in;
}
template<class R> void delete_range(R& r)
{
	for (typename range_iterator<R>::type it = boost::begin(r); it != boost::end(r); ++it)
	{
		delete *it;
	}
}
template<class R1, class R2> static size_t range_min_size (const R1& a, const R2& b)
{
	return min(boost::size(a), boost::size(b));
}
template<class R1, class R2, class R3> static size_t range_min_size (const R1& a, const R2& b, const R3& c)
{
	return min(min(boost::size(a), boost::size(b)), boost::size(c));
}
template<class R1, class R2, class R3, class R4> static size_t range_min_size (const R1& a, const R2& b, const R3& c, const R4& d)
{
	return min(min(min(boost::size(a), boost::size(b)), boost::size(c)), boost::size(d));
}
template<class R1, class R2, class R3, class R4, class R5> static size_t range_min_size (const R1& a, const R2& b, const R3& c, const R4& d, const R5& e)
{
	return min(min(min(min(boost::size(a), boost::size(b)), boost::size(c)), boost::size(d)), boost::size(e));
}
template <class R> static int max_index(const R& r)
{
	return distance(boost::begin(r), max_element(boost::begin(r), boost::end(r)));
}
//ARITHMETIC RANGE OPERATIONS
template<class R1, class R2> static typename range_value<R1>::type inner_product(const R1& a, const R2& b, typename range_value<R1>::type c = 0)
{
	return std::inner_product(boost::begin(a), boost::end(a), boost::begin(b), c);
}
template <class R> static typename range_value<R>::type magnitude(const R& r)
{
	return 0.5 * inner_product(r, r);
}
template <class R1, class R2> static typename range_value<R1>::type sum_of_squares(const R1& r1, const R2& r2)
{
	typename range_const_iterator<R1>::type it1 = boost::begin(r1); 
	typename range_const_iterator<R2>::type it2 = boost::begin(r2); 
	typename range_const_iterator<R1>::type e = boost::end(r1);
	typename range_value<R1>::type v = 0;
	for (; it1 != e; ++it1, ++it2)
	{
		typename range_value<R1>::type diff = *it1 - *it2;
		v += diff * diff;
	}
	return v / 2;
}
template <class R> static typename range_value<R>::type product(const R& r)
{
	return accumulate(boost::begin(r), boost::end(r), (typename range_value<R>::type)1, multiplies<typename range_value<R>::type>());
}
template <class R> static typename range_value<R>::type sum(const R& r)
{
	return accumulate(boost::begin(r), boost::end(r), (typename range_value<R>::type)0);
}
template <class R> static typename range_value<R>::type mean(const R& r)
{
	return sum(r) / (typename range_value<R>::type)boost::size(r);
}
//plus
template<class R1, class R2, class R3> static R1& range_plus(R1& a, const R2& b, const R3& c)
{
	transform(boost::begin(b), boost::end(b), boost::begin(c), boost::begin(a), plus<typename boost::range_value<R1>::type>());
	return a;
}
template<class R1, class R2> static void range_plus_equals(R1& a, const R2& b)
{
	range_plus(a, a, b);
}
//minus
template<class R1, class R2, class R3> static void range_minus(R1& a, const R2& b, const R3& c)
{
	transform(boost::begin(b), boost::end(b), boost::begin(c), boost::begin(a), minus<typename boost::range_value<R1>::type>());
}
template<class R1, class R2> static void range_minus_equals(R1& a, const R2& b)
{
	range_minus(a, a, b);
}
//multiply
template<class R1, class R2> static void range_multiply_val(R1& a, const R2& b, const typename boost::range_value<R2>::type& c)
{
	transform(boost::begin(b), boost::end(b), boost::begin(a), bind2nd(multiplies<typename boost::range_value<R2>::type>(), c));
}
template<class R> static void range_multiply_val(R& a, const typename boost::range_value<R>::type& b)
{
	range_multiply_val(a, a, b);
}
template<class R1, class R2, class R3> static void range_multiply(R1& a, const R2& b, const R3& c)
{
	transform(boost::begin(b), boost::begin(b) + range_min_size(a, b, c), boost::begin(c), boost::begin(a), multiplies<typename boost::range_value<R1>::type>());
}
template<class R1, class R2> static void range_multiply_equals(R1& a, const R2& b)
{
	range_multiply(a, a, b);
}
//divide
template<class R1, class R2> static void range_divide_val(R1& a, const R2& b, const typename boost::range_value<R1>::type& c)
{
	transform(boost::begin(b), boost::end(b), boost::begin(a), bind2nd(divides<typename boost::range_value<R1>::type>(), c));
}
template<class R> static void range_divide_val(R& a, const typename boost::range_value<R>::type& b)
{
	range_divide_val(a, a, b);
}
template<class R1, class R2, class R3> static void range_divide(R1& a, const R2& b, const R3& c)
{
	transform(boost::begin(b), boost::end(b), boost::begin(c), boost::begin(a), divides<typename boost::range_value<R1>::type>());
}
template<class R1, class R2> static void range_divide_equals(R1& a, const R2& b)
{
	range_divide(a, a, b);
}
//SET OPERATIONS
template<class R, class T> void operator +=(set<T>& s, const R& r)
{
	s.insert(boost::begin(r), boost::end(r));
}
//VECTOR OPERATIONS
template<class R, class T> void vector_assign(const R& r, std::vector<T>& v)
{
	v.resize(boost::size(r));
	copy(r, v);
}
//TUPLE OPERATIONS
template<class T1, class T2> static ostream& operator << (ostream& out, const tuple<T1, T2>& t)
{
	out << t.get<0>() << " " << t.get<1>();
	return out;
}
template<class T1, class T2, class T3> static ostream& operator << (ostream& out, const tuple<T1, T2, T3>& t)
{
	out << t.get<0>() << " " << t.get<1>() << " " << t.get<2>();
	return out;
}
template<class T1, class T2, class T3, class T4> static ostream& operator << (ostream& out, const tuple<T1, T2, T3, T4>& t)
{
	out << t.get<0>() << " " << t.get<1>() << " " << t.get<2>() << " " << t.get<3>();
	return out;
}
template<class T1, class T2, class T3, class T4, class T5> static ostream& operator << (ostream& out, const tuple<T1, T2, T3, T4, T5>& t)
{
	out << t.get<0>() << " " << t.get<1>() << " " << t.get<2>() << " " << t.get<3>() << " " << t.get<4>();
	return out;
}
//PAIR OPERATIONS
template<class T1, class T2> static void operator+= (std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
{
	a.first += b.first;
	a.second += b.second;
}
template<class T1, class T2, class T3> static std::pair<T1, T2> operator+ (const std::pair<T1, T2>& a, const T3& b)
{
	return make_pair(a.first + b, a.second + b);
}
template<class T1, class T2> static ostream& operator << (ostream& out, const std::pair<T1, T2>& p)
{
	out << p.first << " " << p.second;
	return out;
}
template<class T1, class T2> static double pair_product(const std::pair<T1, T2>& p)
{
	return (double)(p.first * p.second);
}
template<class T1, class T2> static double pair_sum(const std::pair<T1, T2>& p)
{
	return (double)(p.first + p.second);
}
template<class T1, class T2> static double pair_mean(const std::pair<T1, T2>& p)
{
	return pair_sum(p)/2.0;
}
template <class T1, class T2> static size_t difference(const std::pair<T1,T2>& p)
{
	return p.second - p.first;
}
//MAP OPERATIONS
template<class T1, class T2> static bool in (const std::map<T1, T2>& a, const T1& b)
{
	return (a.find(b) != a.end());
}
template<class T1, class T2> static const T2& at (const std::map<T1, T2>& a, const T1& b)
{
	typename std::map<T1, T2>::const_iterator it = a.find(b);
	check(it != a.end(), str(b) + " not found in map:\n" + str(a));
	return it->second;
}
template<class T1, class T2> static ostream& operator << (ostream& out, const std::map<T1, T2>& m)
{
	for (typename std::map<T1, T2>::const_iterator it = m.begin(); it != m.end(); ++it)
	{
		out << *it << endl; 
	}
	return out;
}
template<class T1, class T2> static ostream& operator << (ostream& out, const std::map<T1, T2*>& m)
{
	for (typename std::map<T1, T2*>::const_iterator it = m.begin(); it != m.end(); ++it)
	{
		out << it->first << " " << *(it->second) << endl; 
	}
	return out;
}
template<class T1, class T2> static T2 sum_right (const std::map<T1, T2>& m)
{
	T2 ret = 0;
	for (typename std::map<T1, T2>::const_iterator it = m.begin(); it != m.end(); ++it)
	{
		ret += it->second;
	}
	return ret;
}
template<class T1, class T2, class T3, class T4> static void operator += (std::map<T1, T2>& a, const std::map<T3, T4>& b)
{
	for (typename std::map<T3, T4>::const_iterator it = b.begin(); it != b.end(); ++it)
	{
		a[it->first] += it->second;
	}
}
template<class T1, class T2, class T3, class T4> static void operator-= (std::map<T1, T2>& a, const std::map<T3, T4>& b)
{
	for (typename std::map<T3, T4>::const_iterator it = b.begin(); it != b.end(); ++it)
	{
		a[it->first] -= it->second;
	}
}
template<class T1, class T2, class T3, class T4> static void operator/= (std::map<T1, T2>& a, const std::map<T3, T4>& b)
{
	for (typename std::map<T3, T4>::const_iterator it = b.begin(); it != b.end(); ++it)
	{
		a[it->first] /= it->second;
	}
}
template<class T1, class T2, class T3, class T4> static void operator*= (std::map<T1, T2>& a, const std::map<T3, T4>& b)
{
	for (typename std::map<T3, T4>::const_iterator it = b.begin(); it != b.end(); ++it)
	{
		a[it->first] *= it->second;
	}
}
template<class T1, class T2, class T3> static void operator*= (std::map<T1, T2>& a, const T3& b)
{
	for (typename std::map<T1, T2>::iterator it = a.begin(); it != a.end(); ++it)
	{
		it->second *= b;
	}
}
template<class T1, class T2, class T3> static void operator/= (std::map<T1, T2>& a, const T3& b)
{
	for (typename std::map<T1, T2>::iterator it = a.begin(); it != a.end(); ++it)
	{
		it->second /= b;
	}
}
template<class R> void delete_map(R& r)
{
	for (typename range_iterator<R>::type it = boost::begin(r); it != boost::end(r); ++it)
	{
		delete it->second;
	}
}
//MULTIMAP OPERATIONS
template<class T1, class T2> static bool in (const multimap<T1, T2>& a, const T1& b)
{
	return (a.find(b) != a.end());
}
//BIMAP OPERATIONS
template<class T1, class T2, class T3, class T4> static ostream& operator << (ostream& out, const boost::bimaps::relation::structured_pair<T1, T2, T3, T4>& p)
{
	out << p.first << " " << p.second;
	return out;
}
template<class T1, class T2> void print_bimap (const bimap<T1, T2>& m, ostream& out)
{
	for (typename bimap<T1, T2>::left_const_iterator it = m.left.begin(); it != m.left.end(); ++it)
	{
		out << *it << endl;
	}
}
template<class T1, class T2> static ostream& operator << (ostream& out, const bimap<T1, T2>& m)
{
	print_bimap(m, out);
	return out;
}
template<class T1, class T2> static bool in_left(const bimap<T1, T2>& a, const T1& b)
{
	return (a.left.find(b) != a.left.end());
}
template<class T1, class T2> static bool in_right(const bimap<T1, T2>& a, const T2& b)
{
	return (a.right.find(b) != a.right.end());
}
//IO OPERATIONS
template<class T> static void print(const T& t, ostream& out = cout)
{
	out << t << endl;
}
template<class T1, class T2> static void print(const T1& t1, const T2& t2, ostream& out = cout)
{
	out << t1 << " " << t2 << endl;
}
template<class T1, class T2, class T3> static void print(const T1& t1, const T2& t2, const T3& t3, ostream& out = cout)
{
	out << t1 << " " << t2 << " " << t3 << endl;
}
template<class T1, class T2, class T3, class T4> static void print(const T1& t1, const T2& t2, const T3& t3, const T4& t4, ostream& out = cout)
{
	out << t1 << " " << t2 << " " << t3  << " " << t4 << endl;
}
template<class T1, class T2, class T3, class T4, class T5> static void print(const T1& t1, const T2& t2, const T3& t3, const T4& t4, const T5& t5, ostream& out = cout)
{
	out << t1 << " " << t2 << " " << t3  << " " << t4 << " " << t5 << endl;
}
static void prt_line(ostream& out = cout)
{
	out << "------------------------------" << endl;
}
template<class T> T read(const std::string& data)
{
	T val;
	std::stringstream ss;
	ss << boolalpha << data;
	check(ss >> val, "cannot read string '" + data + "' into variable with type '" + typeid(T).name() + "'");
	return val;
}
//STRING OPERATIONS
static std::string ordinal(size_t n)
{
	std::string s = str(n);
	if (n < 100)
	{
		char c = nth_last(s);
		if(c == '1')
		{
			return s + "st";
		}
		else if(c == '2')
		{
			return s + "nd";
		}
		else if(c == '3')
		{
			return s + "rd";
		}
	}
	return s + "th";
}
static void trim(std::string& str)
{
    size_t startpos = str.find_first_not_of(" \t\n");
    size_t endpos = str.find_last_not_of(" \t\n");
	if(std::string::npos == startpos || std::string::npos == endpos)  
    {  
        str = "";  
    }  
    else 
	{
        str = str.substr(startpos, endpos-startpos + 1);
	}
}
static bool in(const std::string& str, const std::string& search)
{
	return (str.find(search) != std::string::npos);
}
static bool in(const std::string& str, const char* search)
{
	return in(str, std::string(search));
}
template<class T> std::vector<T>& split(const std::string& original, const char delim = ' ')
{
	static std::vector<T> vect;
	vect.clear();
	std::stringstream ss;
	ss << original;
	std::string s;
	while (getline(ss, s, delim))
	{
		vect += read<T>(s);
	}
	return vect;
}
template<class T, class R >std::string join(const R& r, const std::string joinStr = "")
{
	typename range_iterator<R>::type b = boost::begin(r);
	std::string s = str(*b);
	++b;
	for (; b != end(r); ++b)
	{
		s += joinStr + str(*b);
	}
	return s;
}

//HELPER STRUCTS
template<class T> struct View: public sub_range<std::pair <T*, T*> >
{	
	View(std::pair<T*, T*>& p):
		sub_range<std::pair <T*, T*> >(p)
	{}
	View(T* first = 0, T* second = 0):
		sub_range<std::pair <T*, T*> >(make_pair(first, second))
	{}
	T& at(size_t i)
	{
		check(i < this->size(), "at(" + str(i) + ") called for view of size " + str(this->size()));
		return (*this)[i];
	}
	const T& at(size_t i) const
	{
		check(i < this->size(), "at(" + str(i) + ") called for view of size " + str(this->size()));
		return (*this)[i];	
	}
	template<class T2> const View<T>& operator =(const View<T2>& v) const
	{
		assert(v.size() == this->size());
		copy(v, *this);
		return *this;
	}
};

typedef const tuple<View<double>&, std::vector<double>&>& TVWDVD;
typedef const tuple<View<double>&, View<double>&>& TVWDVWD;

//HACK!! Added by Sergio. For some reason "#include <boost/assign/std/vector.hpp>" does not work
template< class V, class A, class V2 >
inline list_inserter< assign_detail::call_push_back< std::vector<V,A> >, V >
operator+=( std::vector<V,A>& c, V2 v )
{
        return push_back( c )( v );
}



}; /* rnnlib namespace */

#endif
