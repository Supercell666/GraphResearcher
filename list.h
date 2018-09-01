#pragma once
#ifndef _SPRCLIB_LIST
#define _SPRCLIB_LIST

#include <iterator>

namespace sprclib
{
template<class type>
struct list_node
{
    type val;
    struct list_node<type>* left;
    struct list_node<type>* right;
    list_node() :left(nullptr), right(nullptr) {}
    explicit list_node(const type& _val) :val(_val), left(nullptr), right(nullptr) {}
    explicit list_node(type&& _val) :val(std::move(_val)), left(nullptr), right(nullptr) {}
};

template<class type>
class list;

template<class type>
class list_iterator : public std::iterator<std::bidirectional_iterator_tag, type>
{
    friend list<type>;
protected:
    list_node<type>* node;
    list_iterator(list_node<type>* nd) :node(nd) {}
public:
    typedef type& reference;
    typedef type* pointer;
    list_iterator() :node(nullptr) {}
    list_iterator(const list_iterator& it) :node(it.node) {}
    list_iterator(list_iterator&& it) :node(it.node)
    {
        it.node = nullptr;
    }
    list_iterator& operator=(const list_iterator& lst)
    {
        node = lst.node;
        return *this;
    }
    list_iterator& operator=(list_iterator&& lst)
    {
        node = lst.node;
        lst.node = nullptr;
        return *this;
    }
    list_iterator& operator++()
    {
        node = node->right;
        return *this;
    }
    list_iterator& operator++(int)
    {
        auto ret = *this;
        operator++();
        return ret;
    }
    list_iterator& operator--()
    {
        node = node->left;
        return *this;
    }
    list_iterator& operator--(int)
    {
        auto ret = *this;
        operator--();
        return ret;
    }
    reference operator*()
    {
        return node->val;
    }
    const reference operator*() const
    {
        return node->val;
    }
    pointer operator->()
    {
        return &node->val;
    }
    const pointer operator->() const
    {
        return &node->val;
    }
    bool operator==(const list_iterator& it) const
    {
        return node == it.node;
    }
    bool operator!=(const list_iterator& it) const
    {
        return node != it.node;
    }
};

template<class type>
class list_b_iterator :public list_iterator<type>
{
    friend list<type>;
protected:
    list_b_iterator(list_node<type>* nd):list_iterator<type>(nd){}
public:
    list_b_iterator& operator++()
    {
        return list_iterator<type>::operator--();
    }
    list_b_iterator& operator++(int)
    {
        return list_iterator<type>::operator--(1);
    }
    list_b_iterator& operator--()
    {
        return list_iterator<type>::operator++();
    }
    list_b_iterator& operator--(int)
    {
        return list_iterator<type>::operator++(1);
    }
};

template<class type>
class list
{
    friend list_iterator < type >;
private:
    list_node<type> *head, *tail;
    unsigned node_count;
    static const list_iterator<type> _end;
    static const list_b_iterator<type> _rend;

    void _Splice(const list_iterator<type>& it, list& lst)
    {
        auto *before = it.node->left, *after = it.node;
        before->right = lst.head;
        lst.head->left = before;
        lst.tail->right = after;
        after->left = lst.tail;
        node_count += lst.node_count;
        lst.head = lst.tail = nullptr;
        lst.node_count = 0;
    }

    void _Connect(list_node<type> *first_begin, list_node<type> *first_end,
                  list_node<type> *second_begin, list_node<type> *second_end)
    {
        if (first_begin != nullptr)
        {
            first_begin->right = second_begin;
            second_begin->left = first_begin;
        }
        if (first_end != nullptr)
        {
            first_end->left = second_end;
            second_end->right = first_end;
        }
    }
    void _Init(const list& lst)
    {
        if (lst.head == nullptr)
        {
            head = tail = nullptr;
            return;
        }
        head = new list_node<type>();
        head->val = lst.head->val;
        tail = head;
        for (list_node<type>* t = lst.head->right; t != lst.tail->right; t = t->right)
        {
            list_node<type>*node = new list_node<type>();
            tail->right = node;
            tail->right->val = t->val;
            node->left = tail;
            tail = node;
        }
        node_count = lst.node_count;
    }
    void _Init(list&& lst)
    {
        head = lst.head;
        tail = lst.tail;
        node_count = lst.node_count;
        lst.head = nullptr;
        lst.tail = nullptr;
        lst.node_count = 0;
    }

public:
    typedef list_iterator<type> iterator;
    typedef list_b_iterator<type> reverse_iterator;

    list() :head(nullptr), tail(nullptr), node_count(0) {}
    list(const list& lst) : node_count(lst.node_count)
    {
        _Init(lst);
    }
    list(list&& lst)
    {
        _Init(std::move(lst));
    }
    list(const std::initializer_list<type>& lst)
    {
        for(auto& i : lst)
            push_back(i);
    }
    list(std::initializer_list<type>&& lst)
    {
        for(auto&& i : lst)
            push_back(std::move(i));
    }

    void clear()
    {
        while (head != nullptr)
        {
            list_node<type>* t = head->right;
            delete head;
            head = t;
        }
        node_count = 0;
    }
    list& operator=(const list& lst)
    {
        if(this != &lst)
        {
            clear();
            _Init(lst);
        }
        return *this;
    }
    list& operator=(list&& lst)
    {
        if(this != &lst)
        {
            clear();
            _Init(std::move(lst));
        }
        return *this;
    }
    list& operator=(const std::initializer_list<type>& lst)
    {
        clear();
        for(auto& i : lst)
            push_back(i);
        return *this;
    }
    list& operator=(std::initializer_list<type>&& lst)
    {
        clear();
        for(auto&& i : lst)
            push_back(std::move(i));
        return *this;
    }
    iterator erase(const iterator& it)
    {
        if(it.node == tail)
        {
            pop_back();
            return iterator(tail);
        }
        else if(it.node == head)
        {
            pop_front();
            return iterator(head);
        }
        else
        {
            list_node<type>* node = it.node;
            node->left->right = node->right;
            node->right->left = node->left;
            iterator ret(node->left);
            delete node;
            --node_count;
            return ret;
        }
    }
    iterator insert(const iterator& it, const type& val)
    {
        if (it.node == nullptr)
        {
            push_back(val);
            return list_iterator<type>(tail);
        }
        list_node<type> *node = new list_node<type>(std::move(val));
        _Connect(it.node->left, it.node, node, node);
        if (it.node == head)
            head = node;
        if (it.node == tail)
            tail = node;
        ++node_count;
        return list_iterator<type>(node);
    }
    iterator insert(const iterator& it, type&& val)
    {
        if (it.node == nullptr)
        {
            push_back(std::move(val));
            return list_iterator<type>(tail);
        }
        list_node<type> *node = new list_node<type>(std::move(val));
        _Connect(it.node->left, it.node, node, node);
        if (it.node == head)
            head = node;
        if (it.node == tail)
            tail = node;
        ++node_count;
        return list_iterator<type>(node);
    }
    reverse_iterator insert(const reverse_iterator& it, type&& val)
    {
        if (it.node == nullptr)
        {
            push_front(std::move(val));
            return list_b_iterator<type>(tail);
        }
        list_node<type> *node = new list_node<type>(std::move(val));
        _Connect(it.node->left, it.node, node, node);
        if (it.node == head)
            head = node;
        if (it.node == tail)
            tail = node;
        ++node_count;
        return list_b_iterator<type>(node);
    }
    reverse_iterator insert(const reverse_iterator& it, const type& val)
    {
        if (it.node == nullptr)
        {
            push_front(val);
            return list_iterator<type>(tail);
        }
        list_node<type> *node = new list_node<type>(std::move(val));
        _Connect(it.node->left, it.node, node, node);
        if (it.node == head)
            head = node;
        if (it.node == tail)
            tail = node;
        ++node_count;
        return list_iterator<type>(node);
    }
    void push_front(const type& val)
    {
        ++node_count;
        if (head == nullptr)
        {
            head = new list_node<type>(val);
            tail = head;
            return;
        }
        list_node<type> *node = new list_node<type>(val);
        node->right = head;
        head->left = node;
        head = node;
    }
    void pop_front()
    {
        list_node<type> *node = head;
        if (node_count > 1)
            head->right->left = nullptr;
        head = head->right;
        delete node;
        if (--node_count == 0)
            tail = head;
    }
    void push_back(const type& val)
    {
        ++node_count;
        if (head == nullptr)
        {
            head = new list_node<type>(val);
            tail = head;
            return;
        }
        list_node<type> *node = new list_node<type>(val);
        tail->right = node;
        node->left = tail;
        tail = node;
    }

    void pop_back()
    {
        list_node<type> *node = tail;
        if (node_count > 1)
            tail->left->right = nullptr;
        tail = node->left;
        delete node;
        if (--node_count == 0)
            head = tail;
    }
    void splice_back(list& lst)
    {
        if (lst.head == nullptr)
            return;
        if (head == nullptr)
        {
            head = lst.head;
            tail = lst.tail;
            node_count = lst.node_count;
            lst.head = lst.tail = nullptr;
            lst.node_count = 0;
            return;
        }
        tail->right = lst.head;
        lst.head->left = tail;
        tail = lst.tail;
        lst.tail = lst.head = nullptr;
        node_count += lst.node_count;
        lst.node_count = 0;
    }
    void splice_front(list& lst)
    {
        if (lst.head == nullptr)
            return;
        if (tail == nullptr)
        {
            head = lst.head;
            tail = lst.tail;
            node_count = lst.node_count;
            lst.head = lst.tail = nullptr;
            lst.node_count = 0;
            return;
        }
        head->left = lst.tail;
        lst.tail->right = head;
        head = lst.head;
        lst.tail = lst.head = nullptr;
        node_count += lst.node_count;
        lst.node_count = 0;
    }
    void splice(const reverse_iterator& it, list& lst)
    {
        if (it == _rend || it.node == head)
            splice_front(lst);
        else if (it.node == tail)
            splice_back(lst);
        else
            _Splice(it, lst);
    }
    void splice(const iterator& it, list& lst)
    {
        if (it == _end)
            splice_back(lst);
        else if (it.node == head)
            splice_front(lst);
        else
            _Splice(it, lst);
    }
    unsigned size() const
    {
        return node_count;
    }
    iterator begin() const
    {
        return iterator(head);
    }
    const iterator& end() const
    {
        return _end;
    }
    reverse_iterator rbegin() const
    {
        return reverse_iterator(tail);
    }
    const reverse_iterator& rend() const
    {
        return _rend;
    }
    ~list()
    {
        clear();
    }
};
template<class type>
const list_iterator<type> list<type>::_end = nullptr;
template<class type>
const list_b_iterator<type> list<type>::_rend = nullptr;
}

#endif //_SPRCLIB_LIST
