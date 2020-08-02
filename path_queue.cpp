/*************************************
 * Author: M. Babai (M.Babai@rug.nl) *
 * Version:                          *
 * License:                          *
 *************************************/
#include "path_queue.h"

#include <algorithm>
#include <iostream>

PathQueue::PathQueue()
  : m_queueCont(std::vector<int>())
{}

PathQueue::~PathQueue()
{
  m_queueCont.clear();
}

// Copy Constructor
PathQueue::PathQueue (PathQueue const & ot)
  : m_queueCont(ot.m_queueCont)
{}

// Assignment operator.
PathQueue& PathQueue::operator=(PathQueue const &ot)
{
  // Self copy
  if(this != &ot) {
    this->m_queueCont = ot.m_queueCont;
  }
  return (*this);
}

bool PathQueue::isEmpty() const
{
  return (m_queueCont.size() == 0);
}

size_t PathQueue::getNumElement() const
{
  return (this->m_queueCont).size();
}

bool PathQueue::isInQueue(int val) const
{
  std::vector<int>::const_iterator it;
  it = std::find (m_queueCont.begin(), m_queueCont.end(), val);
  return ( it != m_queueCont.end() );
}

void PathQueue::inQueue(int val)
{
  m_queueCont.push_back(val);
}

void PathQueue::deQueue(int val)
{
  std::vector<int>::iterator it;
  it = std::find (m_queueCont.begin(), m_queueCont.end(), val);
  // Empty Queue or element not available
  if(it != m_queueCont.end()) {
    m_queueCont.erase(it);
  }
}

int PathQueue::popFront()
{
  int val = -1;
  std::vector<int>::iterator it = m_queueCont.begin();
  if(m_queueCont.size() != 0) {
    val = *it;
    m_queueCont.erase(it);
  }
  return val;
}

/* Get entry without deleting*/
int PathQueue::GetElement( size_t index) const
{
  if(index < m_queueCont.size()) {
    return (m_queueCont[index]);
  }
  return -1;
}

void PathQueue::Reverse()
{
  std::reverse(m_queueCont.begin(), m_queueCont.end());
}

// Fix this Can go wrong outside the allocated memory
int PathQueue::popBack()
{
  std::vector<int>::iterator it = m_queueCont.end();
  int val = -1;
  if( m_queueCont.size() != 0 ) {
	--it;
    val  = *it;
    m_queueCont.erase(it);
  }
  return val;
}
/* Reset queue. Removes all elements.*/
void PathQueue::ResetQueue()
{
  m_queueCont.clear();
}

void PathQueue::printQueue() const
{
  std::cout << "<INFO> Printing the queue content\n";
  for(size_t i = 0; i < m_queueCont.size(); ++i) {
    std::cout << " " << m_queueCont[i] << ",";
  }
  std::cout << '\n';
}
