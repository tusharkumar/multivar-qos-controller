#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <list>
#include <deque>
#include <set>
#include <string>
#include <ostream>
#include <cassert>
#include <armadillo>
#include <cmath>

namespace fctrl {

	std::string align();
	void incr_align();
	void decr_align();

	double round(double r);

	template<typename T>
	T squared(const T& t)
	{ return t*t; }
	

	template<typename T> class CircularQueue {
	public:
		std::vector<T> q;
		int next;
		int count;

		CircularQueue()
		{
			next = 0;
			count = 0;
		}

		CircularQueue(int length)
		{
			next = 0;
			count = 0;
			q.resize(length);
		}

		void push(const T& t) {
			q.at(next) = t;

			next++;
			if(next == (int)q.size())
				next = 0;
			if(count < (int)q.size())
				count++;
		}

		inline void reset() {
			next = 0;
			count = 0;
		}

		inline void resize_and_reset(int new_length) {
			q.resize(new_length);
			reset();
		}

		inline bool empty() const
		{ return (count != 0); }

		T read_sum() const {
			assert(count > 0);

			T sum = q[0];
			for(int i=1; i<count; i++)
				sum += q[i];

			return sum;
		}

		inline T read_average() const
		{ return read_sum()/count; }

	};


} //namespace fctrl

template<typename T1, typename T2>
	std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& pair)
{
	os << "(" << pair.first << ", " << pair.second << ")";
	return os;
}


template<typename T>
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
	os << "[";
	for(int i=0; i<(int)vec.size(); i++) {
		if(i > 0)
			os << ", ";
		os << vec[i];
	}
	os << "]";

	return os;
}


template<typename T>
	std::ostream& operator<<(std::ostream& os, const std::list<T>& lst)
{
	os << "[";
	for(typename std::list<T>::const_iterator lit = lst.begin(); lit != lst.end(); lit++) {
		if(lit != lst.begin())
			os << " - ";
		os << *lit;
	}
	os << "]";

	return os;
}

template<typename T>
	std::ostream& operator<<(std::ostream& os, const std::deque<T>& vec)
{
	os << "[";
	for(int i=0; i<(int)vec.size(); i++) {
		if(i > 0)
			os << ", ";
		os << vec[i];
	}
	os << "]";

	return os;
}

template<typename T>
	std::ostream& operator<<(std::ostream& os, const std::set<T>& set)
{
	os << "[";
	for(typename std::set<T>::const_iterator sit = set.begin(); sit != set.end(); sit++) {
		if(sit != set.begin())
			os << ", ";
		os << *sit;
	}
	os << "]";

	return os;
}

template<typename T>
	std::ostream& operator<<(std::ostream& os, const std::multiset<T>& set)
{
	os << "[";
	for(typename std::multiset<T>::const_iterator sit = set.begin(); sit != set.end(); sit++) {
		if(sit != set.begin())
			os << ", ";
		os << *sit;
	}
	os << "]";

	return os;
}

template<typename T>
	std::ostream& operator<<(std::ostream& os, const fctrl::CircularQueue<T>& cq)
{
	os << "<";
	int i = cq.next - cq.count;
	if(i < 0)
		i += (int)cq.q.size();
	for(int ic=0; ic<cq.count; ic++) {
		if(i >= (int)cq.q.size())
			i -= (int)cq.q.size();
		if(ic > 0)
			os << ", ";
		os << cq.q.at(i);
		i++;
	}
	os << ">";

	return os;
}

#endif //__UTILS_H__
