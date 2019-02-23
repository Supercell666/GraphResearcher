#pragma once
#ifndef _STRCHANGE_H
#define _STRCHANGE_H

#include <string>
#include <functional>

namespace sprclib
{
	// Вычичсление хэша строки
	long long get_hash_code(const std::string& str)
	{
		long long num1 = 5381;
		long long num2 = num1;
		const char* chPtr2 = str.c_str();
		long long num3, num4;
		while ((num3 = (long long)*chPtr2) != 0)
		{
			num1 = (num1 << 5) + num1 ^ num3;
			num4 = (long long)chPtr2[1];
			if (num4 != 0)
			{
				num2 = (num2 << 5) + num2 ^ num4;
				chPtr2 += 2;
			}
			else
				break;
		}
		return num1 + num2 * 1566083941;
	}

	template<class str_type1, class str_type2>
	// Сравнение подстроки с подстрокой(в алфавитном порядке)
	int substr_comp(const str_type1& str1, unsigned beg1, unsigned end1,
		const str_type2& str2, unsigned beg2, unsigned end2)
	{
		unsigned i, j;
		for (i = beg1, j = beg2; i < end1 && j < end2; ++i, ++j)
		{
			if (str1[i] != str2[j])
				return (str1[i] < str2[j]) ? 1 : -1;
		}
		int k = end1 - beg1 - end2 + beg2;
		if (k < 0)
			return 1;
		else if (k == 0)
			return 0;
		else
			return -1;
	}

	// Типы разделителей
	typedef enum {
		// Между полями строго 1 разделитель
		single,
		// Между полями множество разделителей
		multiplex
	} delimiter_types;

	template<class str_type>
	class pars_params
	{
	public:
		str_type& str;
		bool& need_end;
		unsigned index;

		pars_params(str_type& _str, bool& _need_end, unsigned _index) :
			str(_str), need_end(_need_end), index(_index) {}
		pars_params(const pars_params&) = delete;
		pars_params(pars_params&&) = delete;
	};

#define PARS_ARGS_TYPES(str_type) const sprclib::pars_params<str_type>&
#define PARS_ARGS(str_type) const sprclib::pars_params<str_type>& params
#define SET_PARS_ARGS params
#define BREAKPARS { params.need_end = true; return; }

	// Переводит строку std::string в массив слов std::string
	// str - разбиваемая строка; delim - разделитель;
	// inserter - функция обработки подстрок (принимает pars_args)
	template<class str_type, class inserter_type>
	void pars(const str_type& str, const str_type& delim,
		const inserter_type& inserter,
		const delimiter_types& del_t = delimiter_types::multiplex)
	{
		bool need_end = false; // значение меняется внутри лямбды
		unsigned index = 0; // индекс подстроки(при делении исход. строки разделителем)
		switch (del_t)
		{
		case delimiter_types::multiplex:
			for (int point1 = 0, point2;; point1 = point2 + delim.size())
			{
				while (point1 + delim.size() < str.size() &&
					sprclib::substr_comp(str, point1, point1 + delim.size(),
						delim, 0, delim.size()) == 0)
				{
					point1 += delim.size();
				}
				if (point1 == str.size())
					return;
				point2 = str.find(delim, point1);
				if (point2 == -1)
					point2 = str.size();
				str_type str_set = str.substr(point1, point2 -
					point1);
				inserter(pars_params<str_type>(str_set, need_end, index++));
				if (point2 == str.size() || need_end)
					return;
			}
			break;
		case delimiter_types::single:
			for (int point1 = 0, point2; point1 < str.size();
				point1 = point2 + delim.size())
			{
				point2 = str.find(delim, point1);
				if (point2 == -1)
					point2 = str.size();
				str_type str_set(str.substr(point1, point2 -
					point1));
				inserter(pars_params<str_type>(str_set, need_end, index++));
				if (need_end)
					return;
			}
			break;
		default:
			break;
		}
	}

	// Переводит строку std::string в массив слов std::string
	// str - разбиваемая строка; delim - разделитель;
	// inserter - функция обработки подстрок (принимает pars_args)
	template<class str_type1, class str_type2, class inserter_type>
	inline void pars(const str_type1& str, const str_type2& delim,
		const inserter_type& inserter,
		const delimiter_types& del_t = delimiter_types::multiplex)
	{
		pars(str, str_type1(delim), inserter, del_t);
	}
}
#endif // _STRCHANGE_H
