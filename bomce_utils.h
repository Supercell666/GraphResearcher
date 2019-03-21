#pragma once
#ifndef _BOMCE_UTILS
#define _BOMCE_UTILS

#include <vector>
#include "forward_list.h"
#include "list.h"
#include "bomce.h"

namespace webgr
{
	template<class v_name_type>
    struct cl_vertex;

	template<class sum_type>
	struct bomce_cluster;

	template<class v_name_type>
	using bomce_tmpresult_type = std::vector<sprclib::frwd_list<cl_vertex<
		v_name_type>>>;

	typedef std::vector<float> bomce_params;

	template<class e_count_type>
	using cluster_iterator = typename sprclib::list<
		bomce_cluster<e_count_type>>::iterator;

	template<class sum_type>
	struct cl_cmp
	{
		bool operator()(const cluster_iterator<sum_type>& cl1,
			const cluster_iterator<sum_type>& cl2) const
		{
			return (cl1->id < cl2->id);
		}
	};

	// Паттерн посредник. Без него не возможна "перемотка назад" при перемещении вершин по кластерам
	/*template<class e_count_type>
	struct intermediator
	{
		intermediator<e_count_type>* p; // указатель либо на посредника, либо на итератор на кластер

		intermediator(void* _p = nullptr) :p(_p) {}
		cluster_iterator<e_count_type>& get_it()
		{
			return *dynamic_cast<cluster_iterator<e_count_type>*>(p);
		}
		const cluster_iterator<e_count_type>& get_it() const
		{
			return *dynamic_cast<cluster_iterator<e_count_type>*>(p);
		}
		bomce_cluster<e_count_type>& get()
		{
			return *get_it();
		}
		const bomce_cluster<e_count_type>& get() const
		{
			return *get_it();
		}

		bomce_cluster<e_count_type>& get_next()
		{
			return p->get();
		}
		const bomce_cluster<e_count_type>& get_next() const
		{
			return p->get();
		}

		cluster_iterator<e_count_type>& get_next_it()
		{
			return p->get_it();
		}
		const cluster_iterator<e_count_type>& get_next_it() const
		{
			return p->get_it();
		}

		static intermediator<e_count_type>* new_intermediator(const cluster_iterator<
			e_count_type>& it) // Создаёт новую цепочку ret->p->it и возвращает ret
		{
			return new intermediator<e_count_type>(new intermediator<e_count_type>(
				new cluster_iterator<e_count_type>(it)));
		}
		static intermediator<e_count_type>* new_intermediator(cluster_iterator<
			e_count_type>&& it) // Создаёт новую цепочку ret->p->it и возвращает ret
		{
			return new intermediator<e_count_type>(new intermediator<e_count_type>(
				new cluster_iterator<e_count_type>(std::move(it))));
		}

		// Удаляется из lst по итератору, посредник остаётся
		void erase_next(sprclib::list<bomce_cluster<e_count_type>>& lst)
		{
			if (p->p != nullptr)
			{
				lst.erase(*p->p);
				delete p->p;// указатель на итератор
				p->p = nullptr;
			}			 
		//this->p->it_p
		}
		// Удаляется из lst и разрушает цепочку this->p->it
		void erase(sprclib::list<bomce_cluster<e_count_type>>& lst)
		{
			if (p != nullptr)
			{
				eraze_next();
				delete p; // указатель на посредника
				p = nullptr;
			}
			//this->p->it_p
		}
		// создаём посредника
		void set_next(const cluster_iterator<e_count_type>& it)
		{
			p->p = new cluster_iterator<e_count_type>(it);
		}
		// создаём посредника
		void set_next(cluster_iterator<e_count_type>&& it)
		{
			p->p = new cluster_iterator<e_count_type>(std::move(it));
		}
		bool operator<(const intermediator<e_count_type>& im) const
		{
			return this->get_next()->id < im.get_next()->id;
		}
		bool operator==(const intermediator<e_count_type>& im) const
		{
			return this->get_next()->id == im.get_next()->id;
		}
		bool operator!=(const intermediator<e_count_type>& im) const
		{
			return !(*this == im);
		}
		bool operator>(const intermediator<e_count_type>& im) const
		{
			return im < *this;
		}
		bool operator<=(const intermediator<e_count_type>& im) const
		{
			return !(*this > im);
		}
		bool operator>=(const intermediator<e_count_type>& im) const
		{
			return !(*this < im);
		}
	};

	template<class e_count_type>
	struct v_movement
	{
		intermediator<e_count_type> cl_from, cl_to;
		e_count_type e_count_from, e_count_to;
	};
	*/
}

#endif // _BOMCE_UTILS
