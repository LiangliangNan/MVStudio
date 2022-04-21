#include "match_table.h"


namespace sfm {

		MatchTable::MatchTable(int num_images) {
			match_lists_.resize(num_images);
		}

		void MatchTable::set_match(MatchIndex idx) {
			if (contains(idx))
				return;  // already set

			/* Using vector */
			AdjListElem e;
			e.index = idx.second;
			MatchAdjList &l = match_lists_[idx.first];
			MatchAdjList::iterator p = std::lower_bound(l.begin(), l.end(), e);
			l.insert(p, e);
		}

		void MatchTable::add_match(MatchIndex idx, KeypointMatch m) {
			assert(contains(idx));
			match_list(idx).push_back(m);
		}

		void MatchTable::clear_match(MatchIndex idx) { // but don't erase!
			if (contains(idx)) {
				match_list(idx).clear();
			}
		}

		void MatchTable::remove_match(MatchIndex idx) {
			if (contains(idx)) {
				std::vector<KeypointMatch> &list = match_list(idx);
				list.clear();

				// Remove the neighbor
				AdjListElem e;
				e.index = idx.second;
				MatchAdjList &l = match_lists_[idx.first];
				std::pair<MatchAdjList::iterator, MatchAdjList::iterator> p = 
					std::equal_range(l.begin(), l.end(), e);
				assert(p.first != p.second); // l.end());
				l.erase(p.first, p.second);
			}
		}

		unsigned int MatchTable::num_of_matches(MatchIndex idx) const {
			if (!contains(idx))
				return 0;

			return (unsigned int)match_list(idx).size();
		}

		std::vector<KeypointMatch>& MatchTable::match_list(MatchIndex idx) {
			AdjListElem e;
			e.index = idx.second;
			MatchAdjList &l = match_lists_[idx.first];
			std::pair<MatchAdjList::iterator, MatchAdjList::iterator> p = 
				std::equal_range(l.begin(), l.end(), e);
			assert(p.first != p.second); // l.end());

			return (p.first)->match_list;
		}

		const std::vector<KeypointMatch>& MatchTable::match_list(MatchIndex idx) const {
			AdjListElem e;
			e.index = idx.second;
			const MatchAdjList &l = match_lists_[idx.first];
			std::pair<MatchAdjList::const_iterator, MatchAdjList::const_iterator> p = 
				std::equal_range(l.begin(), l.end(), e);
			assert(p.first != p.second); // l.end());

			return (p.first)->match_list;
		}

		bool MatchTable::contains(MatchIndex idx) const {
			AdjListElem e;
			e.index = idx.second;
			const MatchAdjList &l = match_lists_[idx.first];
			std::pair<MatchAdjList::const_iterator, MatchAdjList::const_iterator> p = 
				std::equal_range(l.begin(), l.end(), e);
			return (p.first != p.second); // l.end());
		}

		void MatchTable::remove_all() {
			for (std::size_t i = 0; i < match_lists_.size(); i++) {
				match_lists_[i].clear();
			}
		}

		unsigned int MatchTable::num_of_neighbors(unsigned int i) {
			return (unsigned int)match_lists_[i].size();
		}

		MatchTable::MatchAdjList& MatchTable::neighbors(unsigned int i) {
			return match_lists_[i];
		}

		MatchTable::MatchAdjList::iterator MatchTable::begin(unsigned int i) {
			return match_lists_[i].begin();
		}

		MatchTable::MatchAdjList::iterator MatchTable::end(unsigned int i) {
			return match_lists_[i].end();
		}

}

