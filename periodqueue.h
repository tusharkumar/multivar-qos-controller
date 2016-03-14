#ifndef __PERIOD_QUEUE_H__
#define __PERIOD_QUEUE_H__

#include <iostream>
#include <queue>
#include <cassert>

namespace fctrl {

	template<typename Data, typename TimeStamp, bool isMinQ = true> class PeriodQueue {
	public:
		class Element {
		public:
			Data data;
			TimeStamp timestamp;

			Element()
			{}

			Element(Data data, TimeStamp timestamp)
				: data(data), timestamp(timestamp)
			{}

			bool operator<(const Element& e) const
			{ return (this->data < e.data) xor isMinQ; }
		};

		std::priority_queue<Element> pq;

		void clear()
		{ pq = std::priority_queue<Element>(); }

		bool empty(TimeStamp cutoff_timestamp) {
			cleantop(cutoff_timestamp);
			return pq.empty();
		}

		Element top(TimeStamp cutoff_timestamp) {
			cleantop(cutoff_timestamp);
			assert(!pq.empty());
			return pq.top();
		}

		inline void push(Element e)
		{ pq.push(e); }

		inline void push(Data data, TimeStamp timestamp)
		{ pq.push( Element(data, timestamp) ); }

		void pop(TimeStamp cutoff_timestamp) {
			cleantop(cutoff_timestamp);
			assert(!pq.empty());
			pq.pop();
		}
	
	private:
		void cleantop(long long int cutoff_timestamp) {
			while(!pq.empty()) {
				Element t = pq.top();
				if(t.timestamp >= cutoff_timestamp)
					break;
				pq.pop();
			}
		}

	};


} //namespace fctrl


template<typename Data, typename TimeStamp, bool isMinQ>
	std::ostream& operator<<(std::ostream& os, const fctrl::PeriodQueue<Data, TimeStamp, isMinQ>& q)
{
	fctrl::PeriodQueue<Data, TimeStamp, isMinQ> copy_q = q;

	os << "[ ";
	while(!copy_q.empty(-1)) {
		typename fctrl::PeriodQueue<Data, TimeStamp, isMinQ>::Element e = copy_q.top(-1);
		os << "(" << e.data << ", " << e.timestamp << ") ";
		copy_q.pop(-1);
	}
	os << "]";
	
	return os;
}



#endif //__PERIOD_QUEUE_H__
