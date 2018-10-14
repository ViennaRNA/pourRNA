
#include "PriorityQueue.h"
#include "StructureUtils.h"

//////////////////   PRIORITY QUEUE KEY  /////////////////////

inline
QueueKey::QueueKey () :
    QueueState ()
{
}

inline
QueueKey::QueueKey (const MyState& state) :
    QueueState (state)
{
}

inline
QueueKey::~QueueKey ()
{
}

inline
QueueValue::QueueValue () :
    QueueState (), StateID (0)
{
}

inline
QueueValue::~QueueValue ()
{
}

inline
bool
QueueKey::operator < (const QueueKey& k2) const
{
  if (QueueState.energy < k2.QueueState.energy)
    return true;
  if (QueueState.energy == k2.QueueState.energy
      && StructureUtils::IsSmaller (QueueState.structure, k2.QueueState.structure)) // use lexicographic order as tiebreaker
    return true;
  return false;
}

inline
bool
QueueKey::operator == (const QueueKey& k2) const
{
  return StructureUtils::IsEqual(QueueState.structure, k2.QueueState.structure); // use lexicographic order as tiebreaker
}

//////////////////   PRIORITY QUEUE IMPLEMENTATION  /////////////////////

template<class QV, class QK, typename CompareKey>
  inline
  PriorityQueue<QV, QK, CompareKey>::PriorityQueue () :
      pq ()
  {
  }

template<class QV, class QK, typename CompareKey>
  inline
  PriorityQueue<QV, QK, CompareKey>::~PriorityQueue ()
  {
  }

///////////////////  ITERATOR ACCESS ///////////////////

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::const_iterator
  PriorityQueue<QV, QK, CompareKey>::begin (void) const
  {
    return pq.begin ();
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::const_iterator
  PriorityQueue<QV, QK, CompareKey>::end (void) const
  {
    return pq.end ();
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::iterator
  PriorityQueue<QV, QK, CompareKey>::begin (void)
  {
    return pq.begin ();
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::iterator
  PriorityQueue<QV, QK, CompareKey>::end (void)
  {
    return pq.end ();
  }

///////////////////  TOP ELEMENT HANDLING  ///////////////////

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::const_iterator
  PriorityQueue<QV, QK, CompareKey>::top (void) const
  {
    return pq.begin (); // first element of the queue according to the order
  }

template<class QV, class QK, typename CompareKey>
  inline
  void
  PriorityQueue<QV, QK, CompareKey>::pop (void)
  {
    // remove first element
    if (pq.begin () != pq.end ())
      {
	pq.erase (pq.begin ());
      }
  }

///////////////////  MISCELLANEOUS FUNCTIONS  ///////////////////

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::size_type
  PriorityQueue<QV, QK, CompareKey>::size (void) const
  {
    // number of queue elements
    return pq.size ();
  }

template<class QV, class QK, typename CompareKey>
  inline
  bool
  PriorityQueue<QV, QK, CompareKey>::empty (void) const
  {
    return pq.empty ();
  }

///////////////////  SEARCH FUNCTIONS  ///////////////////

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::const_iterator
  PriorityQueue<QV, QK, CompareKey>::find (MyState& s) const
  {
    QK k (s);
    return pq.find (k); // first element of the queue according to the order
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::iterator
  PriorityQueue<QV, QK, CompareKey>::find (MyState& s)
  {
    QK k (s);
    return pq.find (k); // first element of the queue according to the order
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::iterator
  PriorityQueue<QV, QK, CompareKey>::find (const QK& key)
  {
    return pq.find (key); // first element of the queue according to the order
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::const_iterator
  PriorityQueue<QV, QK, CompareKey>::find (const QK& key) const
  {
    return pq.find (key); // first element of the queue according to the order
  }

//////////////  INSERT FUNCTIONS  ////////////////////////////

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::InsertResult
  PriorityQueue<QV, QK, CompareKey>::insert (MyState& s)
  {
    QK k (s);
    // k.QueueState = MyState(s);
    // return iterator to newly inserted empty key/value pair to fill
    // auto var = typename InternalPQtype::value_type (k, QV ());
    std::pair<QK, QV> keyValue = std::pair<QK, QV> (k, QV ());
    PriorityQueue<QV, QK, CompareKey>::InsertResult ir = pq.insert (keyValue);
    return ir;
  }

template<class QV, class QK, typename CompareKey>
  inline typename PriorityQueue<QV, QK, CompareKey>::InsertResult
  PriorityQueue<QV, QK, CompareKey>::insert (const QK& key)
  {
    // return iterator to newly inserted empty key/value pair to fill
    QK k (key);
    return pq.insert (typename InternalPQtype::value_type (k, QV ()));
  }

//////////////  ERASE FUNCTIONS  ////////////////////////////

template<class QV, class QK, typename CompareKey>
  inline
  void
  PriorityQueue<QV, QK, CompareKey>::erase (MyState& s)
  {
    QK k (s);
    // erase the corresponding key/value pair
    pq.erase (k);
  }

template<class QV, class QK, typename CompareKey>
  inline
  void
  PriorityQueue<QV, QK, CompareKey>::erase (const QK& key)
  {
    // erase the corresponding key/value pair
    pq.erase (key);
  }

//////////////  ENERGY DEPENDING FUNCTIONS  ////////////////////////////

template<class QV, class QK, typename CompareKey>
  inline
  int
  PriorityQueue<QV, QK, CompareKey>::getMaxE (void) const
  {
    // energy of the last element
    return pq.rbegin ()->first.QueueState.energy;
  }

template<class QV, class QK, typename CompareKey>
  inline
  void
  PriorityQueue<QV, QK, CompareKey>::reduceMaxE (const int maxE)
  {
    // create dummy key of elements to remove
    QK eraseKey;
    eraseKey.QueueState.energy = maxE;
    // delete all elements with key > eraseKey
    pq.erase (pq.upper_bound (eraseKey), pq.end ());
  }
