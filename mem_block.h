//
// Created by bill on 2017-09-29.
//

#ifndef LABB3_MEM_BLOCK_H
#define LABB3_MEM_BLOCK_H

#include <vector>

#define POW 14
#define BLOCKSIZE (1<<POW)
template <class T>
class mem_pool
{
private:
    struct block
    {
        T data[BLOCKSIZE];
    };
    std::vector<block*> blocks;
    int next_free;
public:
    mem_pool();
    void new_block();
    void clear();
    T* new_element();
    T* operator[](int);
};

template <class T>
inline mem_pool<T>::mem_pool()
{
    next_free = 0;
};

template<class T>
inline void mem_pool<T>::new_block()
{
    blocks.push_back(new block);
}

template<class T>
inline T* mem_pool<T>::new_element()
{
  if((next_free&(BLOCKSIZE-1)==0))
  {
      new_block();
  }
    next_free++;
    return &blocks[next_free>>POW]->data[(next_free&(BLOCKSIZE-1))];
};

template<class T>
inline void mem_pool<T>::clear()
{
    for (int i = 0; i < blocks.size(); ++i) {
        delete blocks[i];
    }
}

template<class T>
inline T* mem_pool<T>::operator[](int i)
{
    return &blocks[i>>POW]->data[(i&(BLOCKSIZE-1))];
};


#endif //LABB3_MEM_BLOCK_H
