
#ifndef _SFM_MATCH_TABLE_H_
#define _SFM_MATCH_TABLE_H_

#include "keys.h"

#include <cassert>
#include <algorithm>
#include <string>

#ifndef WIN32
#include <ext/hash_map>
#include <ext/hash_set>
#else
#include <unordered_map>
#include <unordered_set>
#endif



namespace sfm {

	class Transform
	{
	public:
		/* File IO routines */
		void read_from_file(const std::string& file);
		void write_to_file(const std::string& file);

		/* For object movies */
		double	fmatrix[9];
		double	ematrix[9];

		/* For homographies */
		double	H[9];
		double	inlier_ratio;
		int		num_inliers;
	};


	class MatchIndex : public std::pair < unsigned long, unsigned long >
	{
	public:
		MatchIndex(int img1, int img2) : std::pair<unsigned long, unsigned long>(img1, img2) {}
	};



	class MatchTable
	{
	public:
		class AdjListElem {
		public:
			bool operator < (const AdjListElem &other) const { return index < other.index; }
			unsigned int index;
			std::vector<KeypointMatch> match_list;
		};
		typedef std::vector<AdjListElem> MatchAdjList;

	public:
		MatchTable() { }
		MatchTable(int num_images);

		void set_match(MatchIndex idx);
		void add_match(MatchIndex idx, KeypointMatch m);

		void clear_match(MatchIndex idx);

		void remove_match(MatchIndex idx);

		unsigned int num_of_matches(MatchIndex idx) const;

		std::vector<KeypointMatch>& match_list(MatchIndex idx);
		const std::vector<KeypointMatch>& match_list(MatchIndex idx) const;

		bool contains(MatchIndex idx) const;

		void remove_all();

		unsigned int  num_of_neighbors(unsigned int i);
		MatchAdjList& neighbors(unsigned int i);

		MatchAdjList::iterator begin(unsigned int i);
		MatchAdjList::iterator end(unsigned int i);

	private:
		std::vector<MatchAdjList> match_lists_;
	};

}



#ifdef WIN32
namespace stdext {
	template<>
	class hash_compare < sfm::MatchIndex > {
	public:
		static const size_t bucket_size = 4;
		static const size_t min_buckets = 8;
		size_t operator()(const sfm::MatchIndex &__x) const {
			return __x.first * 1529 + __x.second;
		}

		bool operator()(const sfm::MatchIndex &__x1, const sfm::MatchIndex &__x2) const {
			return (__x1.first < __x2.first) || (__x1.first == __x2.first && __x1.second < __x2.second);
		}
	};
}
#else
namespace __gnu_cxx {
	template<>
	struct hash < MatchIndex > {
		size_t
			operator()(MatchIndex __x) const
		{
			return __x.first * 1529 + __x.second;
		}
	};
}
#endif


#ifndef WIN32
typedef __gnu_cxx::hash_set<int>	HashSetInt;
#else
typedef std::unordered_set<int>		HashSetInt;
#endif


namespace sfm {
	class TransformHashMap : public std::unordered_map < sfm::MatchIndex, sfm::Transform, stdext::hash_compare< sfm::MatchIndex > >
	{
	};
}


#endif 
