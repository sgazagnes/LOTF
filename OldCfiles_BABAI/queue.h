#pragma once
#ifndef QUEUE_H
#define QUEUE_H

typedef struct _PrioQueue {
  unsigned long **m_levels;			/* Array of pixels by levels in the queue */
  unsigned long  m_top;				/* Top pixel */
  int m_num_levels;			/* Number of levels in queue */
} PrioQueue;

PrioQueue *create_prio_queue(unsigned long max_rank);
void insert_prio_queue(PrioQueue *q, unsigned long prefix);
void prio_queue_remove(PrioQueue *q);
void free_prio_queue(PrioQueue *q);
unsigned long num_levels(unsigned long size);

#endif
