#ifndef SIMPLEOBJECTPOOL_H
#define SIMPLEOBJECTPOOL_H

#include <vector>
#include <algorithm>
#include <iterator>

namespace icy {
    template<class T>
    class SimplePool;
}

template<class T>
class icy::SimplePool
{

public:
    SimplePool(int initialSize);
    ~SimplePool();
    SimplePool& operator=(SimplePool&) = delete;

    T* take();
    void release(T* p) {available.push(p);}
    void release(std::vector<T*> &vec);
    void releaseAll();
    void printout(); // for testing

private:
    std::vector<T*> available;      // items that are free to use
    std::vector<T*> registry;       // all items of the pool
};

template<class T>
icy::SimplePool<T>::SimplePool(int initialSize)
{
    registry.reserve(initialSize*2);
    available.reserve(initialSize*2);
    for(int i=0;i<initialSize;i++)
    {
        T* obj = new T;
        available.push_back(obj);
        registry.push_back(obj);
    }
}

template<class T>
icy::SimplePool<T>::~SimplePool()
{
    for(auto &x : registry) delete x;
}

template <class T>
T* icy::SimplePool<T>::take()
{
    T *p;
    if(available.size()==0)
    {
        p = new T;
        registry.push_back(p);
    }
    else
    {
        p = available.back();
        available.pop_back();
    }
    return p;
}

template<class T>
void icy::SimplePool<T>::releaseAll()
{
    available.resize(registry.size());
    std::copy(registry.begin(),registry.end(),available.begin());
}

template<class T>
void icy::SimplePool<T>::release(std::vector<T*> &vec)
{
    std::copy(vec.begin(),vec.end(),std::back_inserter(available));
    vec.clear();
}

#endif
