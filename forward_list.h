#pragma once

#ifndef SPRCLIB_FORWARD_LIST
#define SPRCLIB_FORWARD_LIST

#include <forward_list>

namespace sprclib
{
	template<class type>
	class frwd_list :private std::forward_list<type>
	{
	private:
		unsigned el_count;
	public:
		typedef typename std::forward_list<type>::iterator iterator;
		using std::forward_list<type>::begin;
		using std::forward_list<type>::end;
		using std::forward_list<type>::empty;

		frwd_list(frwd_list&& lst) :std::forward_list<type>(std::move(lst)),
			el_count(lst.el_count)
		{
			lst.el_count = 0;
		}
		frwd_list& operator=(frwd_list&& cnt)
		{
			el_count = cnt.el_count;
			std::forward_list<type>::operator=(std::move(cnt));
			//cnt.el_count = 0;
			return *this;
		}	

		unsigned size() const
		{
			return el_count;
		}

		frwd_list() :std::forward_list<type>(), el_count(0) {}

		frwd_list(const std::initializer_list<type>& lst) :
			std::forward_list<type>(lst), el_count(lst.size())
		{}
		template<typename type2>
		frwd_list(const type2& cnt) : std::forward_list<type>(cnt),
			el_count(cnt.size())
		{}

		template<typename... types>
		frwd_list& operator=(const std::initializer_list<types...>& lst)
		{
			std::forward_list<type>::operator=()(lst);
			el_count = lst.size();
			return *this;
		}
		template<typename type2>
		frwd_list& operator=(const type2& cnt)
		{
			std::forward_list<type>::operator=()(cnt);
			el_count = cnt.size();
			return *this;
		}
		void emplace_after(const typename frwd_list::iterator& it, frwd_list& lst)
		{
			if (lst.el_count == 0)
				return;
			if (el_count == 0)
			{
				el_count = lst.el_count;
				std::forward_list<type>::operator=(std::move(lst));
			}
			else
			{
				el_count += lst.size();
				std::forward_list<type>::splice_after(it, lst);
			}
			lst.el_count = 0;
		}
	};
}


#endif // SPRCLIB_FORWARD_LIST