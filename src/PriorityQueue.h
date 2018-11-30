#ifndef PRIORITYQUEUE_H_
#define PRIORITYQUEUE_H_

#include "MyState.h"
#include <map>
#include <set>

//! the key for an element of the queue
class QueueKey
{
public:
MyState QueueState;

//! less comparison based on E and using the lexicographical order
//! of s as tiebreaker
bool
operator<(const QueueKey& toCompare) const;


bool
operator ==(const QueueKey& toCompare) const;


QueueKey ();
QueueKey (const MyState& state);
~QueueKey ();
};

//! the updatable value stored for a key
class QueueValue
{
public:

MyState QueueState;
size_t  StateID;

QueueValue ();
~QueueValue ();
};

//! Implements a priority queue of compressed states that are ordered
//! according to their energy (using lexicographical order on compressed
//! string representation as tiebreaker).
//! The top() element corresponds to the state with the smallest energy
//! in the queue.
template<class QV                           // QueueValue of PriorityQueue_ as template parameter
         , class QK             = QueueKey  // QueueKey with default
         , typename CompareKey  = std::less<QK> >
// comparison operator
class PriorityQueue
{
public:

// access to the template arguments
typedef QV QueueValue;
typedef QK QueueKey;

protected:
//! shortcut for the internal priority queue data structure
typedef typename std::map<QK, QV, CompareKey> InternalPQtype;

//! the internal priority queue data structure
InternalPQtype pq;

public:

//! constant iterator to a <QK, T> pair of the queue
typedef typename InternalPQtype::const_iterator const_iterator;
//! iterator to a <QK,QueueVal> pair of the queue
typedef typename InternalPQtype::iterator iterator;
//! size type of the queue
typedef typename InternalPQtype::size_type size_type;
//! Result type of insertion operation
typedef typename std::pair<iterator, bool> InsertResult;

//! construction
PriorityQueue ();
//! destruction
virtual
~PriorityQueue ();

//! returns a constant iterator to the top() element of the queue
const_iterator
begin(void) const;


//! returns a constant iterator behind the last element of the queue
const_iterator
end(void) const;


//! returns an iterator to the top() element of the queue
iterator
begin(void);


//! returns an iterator behind the last element of the queue
iterator
end(void);


//! returns the top element of the queue which is the top() element
const_iterator
top(void) const;


//! deletes the top element of the queue
void
pop(void);


//! the number of elements in the queue
//! @return the number of elements
size_type
size(void) const;


//! Whether or not the queue contains no elements
//! @return true if size()==0 and false otherwise
bool
empty(void) const;


//! Returns a constant iterator to the queue element representing s
//! or end() otherwise.
//! @param s the state of which the representation has to be found
//! @return constant iterator to the element or end() if not found
const_iterator
find(MyState& s) const;


//! Returns an iterator to the queue element representing s
//! or end() otherwise.
//! @param s the state of which the representation has to be found
//! @return iterator to the element or end() if not found
iterator
find(MyState& s);


//! Returns a constant iterator to the queue element representing s
//! or end() otherwise.
//! @param key the key to find
//! @return constant iterator to the element or end() if not found
const_iterator
find(const QK& key) const;


//! Returns an iterator to the queue element representing s
//! or end() otherwise.
//! @param key the key to find
//! @return iterator to the element or end() if not found
iterator
find(const QK& key);


//! Inserts a key/value pair to the queue using the given information.
//! A new QueueValue object is created and can be filled afterwards
//! using the provided iterator.
//! NOTE: An already existing QueueValue object for the given key
//! will be overwritten!
//! @param keyState the state that is to be inserted as a key
//! @return (first) the iterator to the (newly created) QueueValue
//!         object and (second) a boolean that is true if the value
//!         was created and NOT already part of the queue.
InsertResult
insert(MyState& keyState);


//! Inserts an empty QueueValue for the given key.
//! The new QueueValue object is created and can be filled afterwards
//! using the provided iterator.
//! NOTE: An already existing QueueValue object for the given key
//! will be overwritten!
//! @param key the key a new value is to be inserted
//! @return (first) the iterator to the (newly created) QueueValue
//!         object and (second) a boolean that is true if the value
//!         was created and NOT already part of the queue.
InsertResult
insert(const QK& key);


//! Erases the key/value pair representing the given state
//! @param s the state of that the key/value pair is to erase
void
erase(MyState& s);


//! Erases the key/value pair of the given key
//! @param key the key to erase the key/value pair of
void
erase(const QK& key);


//! Returns the maximal energy stored in the queue
//! @return the energy of the last queue element
int
getMaxE(void) const;


//! Deletes all elements from the queue with key.E > maxE
//! @param maxE the energy to that all elements of the queue should
//!             be smaller
void
reduceMaxE(const int maxE);
};

#include "PriorityQueue.cpp"

#endif /*PRIORITYQUEUE_H_*/
