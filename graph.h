#pragma once
#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <map>
#include <ostream>
#include <functional>
#include <memory>
#include <queue>
#include <iterator>
#include <type_traits>
#include "graphutility.h"

namespace webgr
{
	template<typename ptr_type>
	class hash_ptr : public std::hash<typename std::remove_pointer<ptr_type>::type>
	{
		typedef typename std::hash<typename std::remove_pointer<ptr_type>::type> parrent;
	public:
		typedef typename parrent::result_type result_type;
		typedef ptr_type argument_type;

		size_t operator()(const argument_type& ptr) const
		{
			return parrent::operator()(*ptr);
		}
	};

	template<typename ptr_type>
	struct compare_ptr
	{
		bool operator()(const ptr_type& first, const ptr_type& second) const
		{
			return (*first < *second);
		}
	};

	template<typename ptr_type>
	struct eq_ptr
	{
		bool operator()(const ptr_type& first, const ptr_type& second) const
		{
			return (*first == *second);
		}
	};

	template<class type, class edge_value_type, class e_count_type = unsigned,
		template<class> class compare_type = compare_ptr>
	class graph;

	template<class type, class edge_value_type = unsigned>
	class vertex
	{
		friend graph<type, edge_value_type>;
	private:
		// Имя вершины
		type* _name;
	public:
		// Входящая степень
		edge_value_type in_d,
			// Исходящая степень
			out_d,
			// Количество петель
			loop;
		// Входящие рёбра
		std::map<unsigned, edge_value_type> input,
			// Исходящие рёбра
			output;

		vertex(type* _Name = nullptr) : in_d(0), out_d(0), loop(0), _name(_Name) {}
		vertex(const vertex& lnk) :input(lnk.input),
			output(lnk.output), in_d(lnk.in_d), out_d(lnk.out_d),
			loop(lnk.loop), _name(new type(*lnk._name)) {}
		vertex(vertex&& lnk) :input(std::move(lnk.input)),
			in_d(std::move(lnk.in_d)), out_d(std::move(lnk.out_d)),
			output(std::move(lnk.output)), loop(std::move(lnk.loop)),
			_name(std::move(lnk._name))
		{
			lnk._name = nullptr;
		}

		vertex& operator=(vertex&& v)
		{
			if (this == &v)
				return *this;
			_name = v._name;
			v._name = nullptr;
			in_d = v.in_d;
			out_d = v.out_d;
			loop = v.loop;
			input = std::move(v.input);
			output = std::move(v.output);
			return *this;
		}
		vertex& operator=(const vertex& v)
		{
			if (this == &v)
				return *this;
			_name = new type(*v._name);
			in_d = v.in_d;
			out_d = v.out_d;
			loop = v.loop;
			input = v.input;
			output = v.output;
			return *this;
		}
		const type& name() const
		{
			return *_name;
		}

		edge_value_type deg() const
		{
			return in_d + out_d;
		}

		template<class _Function>
		void foreach_neighbors(const _Function& func) const
		{
			typename std::map<unsigned, edge_value_type>::const_iterator ob = output.begin(),
				oe = output.end(), ib = input.begin(), ie = input.end();
			merging(ob, oe, ib, ie, func);
		}

		~vertex()
		{
			if (_name != nullptr)
				delete _name;
		}
	};

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type, class iterator_type>
	class graph_iterator : public std::iterator < std::bidirectional_iterator_tag,
		typename std::pair<unsigned const&, vertex<type, edge_value_type>& >>
	{
		friend graph<type, edge_value_type, e_count_type, compare_type>;
	public:
		typedef typename std::pair<unsigned const&,
			vertex<type, edge_value_type> const& > value_type;
		typedef typename std::pair<unsigned const&,
			vertex<type, edge_value_type> const& >* pointer;
		typedef typename std::pair<unsigned const&,
			vertex<type, edge_value_type> const& >& reference;
	protected:
		const graph<type, edge_value_type, e_count_type, compare_type>* gr;
		iterator_type it;
		mutable pointer p; // пара из двух ссылок (lasy паттерн)

		graph_iterator(const graph<type, edge_value_type,
			e_count_type, compare_type>* _gr, const iterator_type& _it) :
			gr(_gr), it(_it), p(nullptr)
		{}
	public:
		graph_iterator() :gr(nullptr), p(nullptr) {}
		graph_iterator(const graph_iterator& iter) :gr(iter.gr), it(iter.it),
			p((iter.p != nullptr) ? (new value_type(*iter.p)) : nullptr) {}
		graph_iterator(graph_iterator&& iter) :gr(iter.gr), it(std::move(iter.it)),
			p(iter.p)
		{
			iter.p = nullptr;
			iter.gr = nullptr;
		}

		graph_iterator& operator=(const graph_iterator& iter)
		{
			if (this == &iter)
				return *this;
			gr = iter.gr;
			it = iter.it;
			p = (iter.p != nullptr) ? (new value_type(*iter.p)) : nullptr;
			return *this;
		}

		graph_iterator& operator=(graph_iterator&& iter)
		{
			if (this == &iter)
				return *this;
			gr = iter.gr;
			it = std::move(iter.it);
			p = iter.p;
			iter.gr = nullptr;
			if (iter.p != nullptr)
			{
				delete iter.p;
				iter.p = nullptr;
			}
			return *this;
		}

		graph_iterator& operator++()
		{
			++it;
			if (p != nullptr)
			{
				delete p;
				p = nullptr;
			}
			return *this;
		}

		graph_iterator operator++(int)
		{
			auto temp = *this;
			this->operator++();
			return temp;
		}

		graph_iterator& operator--()
		{
			--it;
			if (p != nullptr)
			{
				delete p;
				p = nullptr;
			}
			return *this;
		}

		graph_iterator operator--(int)
		{
			auto temp = *this;
			this->operator--();
			return temp;
		}

		template<class gr_iterator>
		bool operator==(const graph_iterator<type, edge_value_type,
			e_count_type, compare_type, gr_iterator> &it) const
		{
			return (it.it == this->it);
		}

		template<class gr_iterator>
		bool operator!=(const graph_iterator<type, edge_value_type,
			e_count_type, compare_type, gr_iterator>& _it) const
		{
			return !this->operator==(_it);
		}

		reference operator*()
		{
			if (p == nullptr)
				p = new value_type(it->second, gr->operator[](it->second));
			return *p;
		}

		const reference operator*() const
		{
			if (p == nullptr)
				p = new value_type(it->second, gr->operator[](it->second));
			return *p;
		}

		pointer operator->()
		{
			if (p == nullptr)
				p = new value_type(it->second, gr->operator[](it->second));
			return p;
		}

		const pointer operator->() const
		{
			if (p == nullptr)
				p = new value_type(it->second, gr->operator[](it->second));
			return p;
		}

		~graph_iterator()
		{
			if (p != nullptr)
			{
				delete p;
				p = nullptr;
			}
		}
	};

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	class graph
	{
		friend graph_iterator<type, edge_value_type, e_count_type, compare_type, typename
			std::map<type*, unsigned, compare_type<type*>>::iterator>;
	public:
		typedef graph_iterator<type, edge_value_type, e_count_type, compare_type, typename
			std::map<type*, unsigned, compare_type<type*>>::iterator> iterator;
		typedef edge_value_type edge_type;
		typedef type vertex_name_type;
		typedef e_count_type edge_count_type;
	private:
		std::queue<unsigned> deleted;
		edge_count_type m; // Число рёбер
		mutable std::vector<vertex<type, edge_value_type>> v; // индекс вершины на её параметры
		mutable std::map<type*, unsigned, compare_type<type*>> reindexed; // имя вершины на её индекс
		mutable iterator* it_end; // не меняется при вставке/удалении

		void erase(type*, unsigned);
	public:
		graph() :m(0), it_end(nullptr) {}
		graph(const graph& gr) :v(gr.v), m(gr.m), reindexed(gr.reindexed),
			deleted(gr.deleted), it_end(nullptr) {}
		graph(graph&& gr) :v(std::move(gr.v)), m(std::move(gr.m)),
			reindexed(std::move(gr.reindexed)), deleted(std::move(gr.deleted))
		{
			gr.m = 0;
			it_end = nullptr;
			gr.it_end = nullptr;
		}

		void operator=(graph&& gr)
		{
			m = gr.m;
			v = std::move(gr.v);
			reindexed = std::move(gr.reindexed);
			it_end = nullptr;
			gr.it_end = nullptr;
			deleted = std::move(gr.deleted);
			gr.m = 0;
		}

		graph& operator=(const graph& gr)
		{
			m = gr.m;
			v = gr.v;
			reindexed = gr.reindexed;
			it_end = nullptr;
			deleted = gr.deleted;
			return *this;
		}

		// добавление вершины с заданным именем
		iterator new_vertex(const type&);
		iterator _New_Vertex(const type&, unsigned);

		// Поиск вешины по имени
		iterator find(const type& v_name) const
		{
			return iterator(this, reindexed.find((type*)&v_name));
		}

		// Удаление вершины
		void erase(const iterator& it)
		{
			erase(*it.it);
		}

		// Удаление вершины
		void erase_by_name(const type& name)
		{
			auto it = this->find(name);
			if (it != this->end())
				erase(*it.it);
		}

		void erase_by_index(unsigned index)
		{
			erase(v[index]._name, index);
		}

		// Установка n связей first --> second
		void connect_by_index(unsigned, unsigned, edge_value_type n = 1);

		// Установка n связей first --> second
		void connect(const iterator& it_f,
			const iterator& it_s, edge_value_type n = 1)
		{
			connect_by_index(it_f.it->second, it_s.it->second, n);
		}

		// Установка n связей first --> second
		void connect_by_name(const type& first, const type& second,
			edge_value_type n = 1)
		{
			auto it_f = reindexed.find((type*)&first);
			if (first == second)
				this->connect_by_index(it_f->second, it_f->second, n);
			else
				this->connect_by_index(it_f->second, reindexed[(type*)&second], n);
		}

		// Удаление n рёбер first --> second
		// Если рёбер меньше, чем n - они удалятся все
		// Возвращает количество удалённых рёбер
		edge_value_type disconnect_by_index(unsigned,
			unsigned, edge_value_type);

		// Удаляет се ребра first --> second
		edge_value_type disconnect_by_index(unsigned,
			unsigned);

		// Удаление n рёбер first --> second
		// Если рёбер меньше, чем n - они удалятся все
		// Возвращает количество удалённых рёбер
		edge_value_type disconnect(const iterator& it_f,
			const iterator& it_s, edge_value_type n)
		{
			return this->disconnect_by_index(it_f.it->second,
				it_s.it->second, n);
		}

		// Удаляет все рёбра it_f --> it_s
		edge_value_type disconnect(const iterator& it_f,
			const iterator& it_s)
		{
			return this->disconnect_by_index(it_f.it->second,
				it_s.it->second);
		}

		// Удаление n рёбер first --> second
		// Если рёбер меньше, чем n - они удалятся все
		// Возвращает количество удалённых рёбер
		edge_value_type disconnect_by_name(const type&,
			const type&, edge_value_type n = 1);

		// Удаляет все ребра first --> second
		edge_value_type disconnect_by_name(const type&,
			const type&);

		// Количество рёбер it_f --> it_s
		edge_value_type count_by_index(unsigned it_f,
			unsigned it_s) const
		{
			if (it_f == it_s)
				return v[it_f].loop;
			auto it = v[it_f].output.find(it_s);
			return (it == v[it_f].output.end()) ? 0 : it->second;
		}

		// Количество рёбер it_f --> it_s
		edge_value_type count(const iterator& it_f,
			const iterator& it_s) const
		{
			return this->count_by_index(it_f.it->second, it_s.it->second);
		}

		// Количество рёбер first --> second
		edge_value_type count_by_name(const type&, const type&) const;

		void merge_by_index(unsigned, unsigned);
		// Объединение двух вершин. it_f += it_s
		void merge(const iterator& it_f, const iterator& it_s)
		{
			merge_by_index(it_f->it->first, it_s->it->first);
		}

		// Объединение двух вершин. first += second
		void merge_by_name(const type& first, const type& second)
		{
			merge(this->find(first), this->find(second));
		}

		// Количество вершин в графе
		unsigned size() const
		{
			return reindexed.size();
		}

		// Физически хранимое число вершин(не изменяется при удалении вершины)
		unsigned capacity() const
		{
			return v.size();
		}

		const vertex<type, edge_value_type>& operator[](unsigned index) const
		{
			return v[index];
		}

		void clear()
		{
			v.clear();
			deleted = std::move(std::queue<unsigned>());
			for (auto& i : reindexed)
				delete i.first;
			reindexed.clear();
			m = 0;
		}

		bool _Resize(size_t _size)
		{
			if (_size < v.size())
				return false;
			v.resize(_size);
			return true;
		}

		// переиндексация вершин графа
		// устраняет незадействованные индексы
		void reindex();

		iterator begin()
		{
			return iterator(this, this->reindexed.begin());
		}

		const iterator begin() const
		{
			return iterator(this, this->reindexed.begin());
		}

		const iterator& end() const
		{
			if (it_end == nullptr)
			{
				it_end = new iterator(this, this->reindexed.end());
			}
			return *it_end;
		}

		// Число рёбер
		edge_count_type edge_count() const
		{
			return m;
		}

		edge_value_type erase_loops();
	};


	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Удаляет петли
	edge_value_type graph<type, edge_value_type, e_count_type,
		compare_type>::erase_loops()
	{
		edge_value_type ret;
		if (std::is_arithmetic<edge_value_type>::value)
			ret = 0;
		for (auto& i : this->v)
		{
			m -= i.loop;
			ret += i.loop;
			i.loop = 0;
		}
		return ret;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// добавление вершины с заданным именем
	typename graph<type, edge_value_type, e_count_type, compare_type>::
		iterator graph<type, edge_value_type, e_count_type,
		compare_type>::new_vertex(const type& name)
	{
		auto it = reindexed.find((type*)&name);
		if (it == reindexed.end())
		{
			auto* p = new type(name);
			if (this->deleted.empty())
			{
				it = reindexed.insert({ p, v.size() }).first;
				v.push_back(vertex<type, edge_value_type>(p));
			}
			else
			{
				it = reindexed.insert({ p, deleted.front() }).first;
				if (v[deleted.front()]._name != nullptr)
					delete v[deleted.front()]._name;
				v[deleted.front()]._name = p;
				deleted.pop();
			}
			if (it_end != nullptr)
			{
				delete it_end;
				it_end = nullptr;
			}
		}
		return graph<type, edge_value_type, e_count_type,
			compare_type>::iterator(this, it);
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	typename graph<type, edge_value_type, e_count_type, compare_type>::
		iterator graph<type, edge_value_type, e_count_type,
		compare_type>::_New_Vertex(const type& name, unsigned index)
	{
		auto it = reindexed.find((type*)&name);
		if (it == reindexed.end())
		{
			auto* p = new type(name);
			it = reindexed.insert({ p, index }).first;
			if (v[index]._name != nullptr)
				delete v[index]._name;
			v[index]._name = p;
			if (it_end != nullptr)
			{
				delete it_end;
				it_end = nullptr;
			}
		}
		return graph<type, edge_value_type, e_count_type,
			compare_type>::iterator(this, it);
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	void graph<type, edge_value_type, e_count_type, compare_type>::erase(
		type* p, unsigned index)
	{
		for (auto& i : v[index].input)
			this->disconnect_by_index(index, i.first, i.second);
		for (auto& i : v[index].output)
			this->disconnect_by_index(i.first, index, i.second);
		if (v[index]._name != nullptr)
		{
			delete v[index]._name;
			v[index]._name = nullptr;
		}
		v[index].in_d = v[index].out_d = v[index].loop = 0;
		m -= v[index].loop;
		deleted.push(index);
		reindexed.erase(p);
		if (it_end != nullptr)
		{
			delete it_end;
			it_end = nullptr;
		}
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Установка n рёбер it_f --> it_s
	void graph<type, edge_value_type, e_count_type, compare_type>::
		connect_by_index(unsigned it_f, unsigned it_s, edge_value_type n)
	{
		if (it_f == it_s)
		{
			this->v[it_f].loop += n;
		}
		else
		{
			auto it = v[it_f].output.find(it_s);
			if (it == v[it_f].output.end())
			{
				v[it_f].output.insert({ it_s, n });
				v[it_s].input.insert({ it_f, n });
			}
			else
			{
				it->second += n;
				v[it_s].input[it_f] += n;
			}
		}
		v[it_f].out_d += n;
		v[it_s].in_d += n;
		m += n;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	edge_value_type graph<type, edge_value_type, e_count_type,
		compare_type>::count_by_name(const type& first,
			const type& second) const
	{
		auto it_f = this->reindexed.find(&first);
		if (it_f == this->reindexed.end())
			return 0;
		if (first == second)
			return this->count_by_index(it_f->second, it_f->second);
		auto it_s = this->reindexed.find(&second);
		if (it_s != this->end())
			return this->count_by_index(it_f->second, it_s->second);
		return 0;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Удаление n рёбер it_f --> it_s
	// Если рёбер меньше, чем n - они удалятся все
	// Возвращает количество удалённых рёбер
	edge_value_type graph<type, edge_value_type, e_count_type,
		compare_type>::disconnect_by_index(unsigned it_f,
			unsigned it_s, edge_value_type n)
	{
		edge_value_type ret = n;
		if (it_f == it_s)
		{
			if (v[it_f].loop > n)
				v[it_f].loop -= n;
			else
			{
				ret = v[it_f].loop;
				v[it_f].loop = 0;
			}
		}
		else
		{
			auto it_fs = v[it_f].output.find(it_s);
			if (it_fs == v[it_f].output.end())
				return 0;
			auto it_sf = v[it_s].input.find(it_f);
			if (it_fs->second > n)
			{
				it_fs->second -= n;
				it_sf->second -= n;
			}
			else
			{
				ret = it_fs->second;
				v[it_s].input.erase(it_sf);
				v[it_f].output.erase(it_fs);
			}
		}
		v[it_f].out_d -= ret;
		v[it_s].in_d -= ret;
		m -= ret;
		return ret;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Удаление n рёбер it_f --> it_s
	// Если рёбер меньше, чем n - они удалятся все
	// Возвращает количество удалённых рёбер
	edge_value_type graph<type, edge_value_type, e_count_type,
		compare_type>::disconnect_by_index(unsigned it_f, unsigned it_s)
	{
		edge_value_type ret;
		if (it_f == it_s)
		{
			ret = v[it_f].loop;
			v[it_f].loop = 0;
		}
		else
		{
			auto it_fs = v[it_f].output.find(it_s);
			if (it_fs == v[it_f].output.end())
				return 0;
			ret = it_fs->second;
			v[it_s].input.erase(it_f);
			v[it_f].output.erase(it_fs);
		}
		v[it_f].out_d -= ret;
		v[it_s].in_d -= ret;
		m -= ret;
		return ret;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Удаление n рёбер first --> second
	// Если рёбер меньше, чем n - они удалятся все
	// Возвращает количество удалённых рёбер
	edge_value_type graph<type, edge_value_type, e_count_type,
		compare_type>::disconnect_by_name(const type& first,
			const type& second, edge_value_type n)
	{
		auto it_f = this->reindexed.find(&first);
		if (it_f != this->reindexed.end())
		{
			if (first == second)
			{
				return this->disconnect_by_index(it_f->second, it_f->second, n);
			}
			else
			{
				auto it_s = this->reindexed.find(&second);
				if (it_s != this->reindexed.end())
				{
					return this->disconnect_by_index(it_f->second,
						it_s->second, n);
				}
			}
		}
		return 0;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Удаление n рёбер first --> second
	// Если рёбер меньше, чем n - они удалятся все
	// Возвращает количество удалённых рёбер
	edge_value_type graph<type, edge_value_type, e_count_type,
		compare_type>::disconnect_by_name(const type& first,
			const type& second)
	{
		auto it_f = this->reindexed.find(&first);
		if (it_f != this->reindexed.end())
		{
			if (first == second)
			{
				return this->disconnect_by_index(it_f->second, it_f->second);
			}
			else
			{
				auto it_s = this->reindexed.find(&second);
				if (it_s != this->reindexed.end())
				{
					return this->disconnect_by_index(it_f->second,
						it_s->second);
				}
			}
		}
		return 0;
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// Объединение двух вершин(second перехоит в first)
	void graph<type, edge_value_type, e_count_type, compare_type>::merge_by_index(
		unsigned first, unsigned second)
	{
		edge_value_type count = this->count_by_index(first, second);
		if (count > 0)
		{
			disconnect_by_index(first, second, count);
			connect_by_index(first, first, count);
		}
		count = this->count_by_index(second, first);
		if (count > 0)
		{
			disconnect_by_index(second, first, count);
			connect_by_index(first, first, count);
		}
		count = this->count_by_index(second, second);
		if (count > 0)
		{
			disconnect_by_index(second, second, count);
			connect_by_index(first, first, count);
		}
		// пересадка рёбер
		for (auto& i : v[second].input)
			connect_by_index(i.first, first, i.second);
		for (auto& i : v[second].output)
			connect_by_index(first, i.first, i.second);
		this->erase_by_index(second);
	}

	template<class type, class edge_value_type, class e_count_type,
		template<class> class compare_type>
	// переиндексация вершин графа
	// устраняет незадействованные индексы
	void graph<type, edge_value_type, e_count_type, compare_type>::reindex()
	{
		unsigned left = 0, right = v.size() - 1;
		for (; left < right; ++left, --right)
		{
			for (; v[left]._name != nullptr; ++left)
			{
				if (left == right)
					goto end;
			}
			for (; v[right]._name == nullptr; --right)
			{
				if (left == right)
					goto end;
			}
			if (left == right)
				goto end;
			v[left] = std::move(v[right]);
			for (auto& i : v[left].input)
			{
				auto t = std::move(v[i].output[right]);
				v[i].output.erase(right);
				v[i].output.insert({ left, std::move(t) });
			}
			for (auto& i : v[left].output)
			{
				auto t = std::move(v[i].input[right]);
				v[i].input.erase(right);
				v[i].input.insert({ left, std::move(t) });
			}
			reindexed[v[left]._name] = left;
		}
	end:
		v.resize(reindexed.size());
		deleted = std::move(std::queue<unsigned>());
	}

	typedef graph<unsigned, unsigned, unsigned, compare_ptr> default_graph;

}

#endif // _GRAPH_H_
