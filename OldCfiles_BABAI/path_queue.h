/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
/*
 * NOTE: This is not very efficient. It is maybe better to implement
 * queue using std::list as insertion and deletion has constant time
 * and goes via link remove and insertion per node. But our data size
 * is very small????
 */
/*
 * Simple Queue (FIFO) implementation. It is not very efficient yet.
 */
#pragma once
#ifndef PATH_QUEUE_CLASS_H
#define PATH_QUEUE_CLASS_H

#include <cstdlib>
#include <vector>

struct PathQueue {
public:
  // Constructor
  PathQueue();

  // Destructor
  virtual ~PathQueue();
  
  // Copy
  PathQueue (PathQueue const & ot);

  // Assign
  PathQueue& operator=(PathQueue const &ot);

  //_____________ Member functions and modifier __________
  /* If empty */
  bool isEmpty() const;  

  /* Number of elements. */
  size_t getNumElement() const;

  /* If an element is in queue*/
  bool isInQueue(int val) const;

  /* Insert element (tail) */
  void inQueue(int val);

  /* Delete element with the given value*/
  void deQueue(int val);

  /* First in the row */
  int  popFront();

  /* Get element without deleting */
  int GetElement( size_t index) const;

  /* Last in the row */
  int  popBack();

  /* Reverse the element order*/
  void Reverse();

  /* Reset queue. Removes all elements.*/
  void ResetQueue();
  
  inline std::vector<int> const& GetListOfElements() const;

  void printQueue() const;

  //________________________ data Members __________
  
  // Container to hold elements(use linked list, is faster for random
  // access and deleting)
  std::vector<int> m_queueCont;

  // protected:
  // private:
};
//_____________________________________________________________
// Function implementations
std::vector<int> const& PathQueue::GetListOfElements() const
{
  return this->m_queueCont;
}
//_____________________________________________________________
#endif// End of interface
