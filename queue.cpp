#include <climits>
#include "queue.h"



int bit_scan_reverse(unsigned long val) {
  return sizeof(val) * CHAR_BIT - 1 - __builtin_clzl(val);
}
int bits_per_word_log2(void) {
  return sizeof(unsigned) * CHAR_BIT - 1 - __builtin_clz((sizeof(unsigned long) * CHAR_BIT));
}
int bits_per_word(void) { 
  return sizeof(unsigned long) * CHAR_BIT;
}

unsigned long num_levels(unsigned long size){
  int bits_per_rank = bit_scan_reverse(size) + 1;
  // number of required levels
  return (bits_per_rank + bits_per_word_log2() - 1) / bits_per_word_log2();
}

PrioQueue *create_prio_queue(unsigned long max_rank) {
  PrioQueue *queue = (PrioQueue *) malloc(sizeof(PrioQueue));
  queue->m_num_levels = num_levels(max_rank);
  queue->m_levels = (unsigned long **) malloc(queue->m_num_levels * sizeof(unsigned long*));
  for (int i = queue->m_num_levels; i--;){
    max_rank >>= bits_per_word_log2();
    queue->m_levels[i] = (unsigned long *) malloc((max_rank+1) * sizeof(unsigned long));
    memset(queue->m_levels[i], 0,  (max_rank+1)*sizeof(unsigned long));
  }
  queue->m_top = 0;
  return queue;
}


void insert_prio_queue(PrioQueue *q, unsigned long prefix) {

  if (prefix > q->m_top) {
    q->m_top = prefix;

  }
   for (int i = q->m_num_levels; i--;) {     
     unsigned long word_idx = prefix >> bits_per_word_log2();
     int bit_idx = prefix & (bits_per_word() - 1);
     q->m_levels[i][word_idx] |= (ulong)(1) << bit_idx;
     prefix = word_idx;
   }
}                              


void prio_queue_remove(PrioQueue *q) {
  unsigned long prefix = q->m_top;
  for (int i = q->m_num_levels; i--;){
    unsigned long word_idx = prefix >> bits_per_word_log2();
    q->m_levels[i][word_idx] &= ~((unsigned long)(1) << (prefix & (bits_per_word() - 1)));
    if (q->m_levels[i][word_idx] != 0) {
      for (int j = i; j != q->m_num_levels; ++j)	
	  word_idx = (word_idx << bits_per_word_log2()) ^ bit_scan_reverse(q->m_levels[j][word_idx]);  
	
      q->m_top = word_idx;

      return;
    }
    prefix = word_idx;
  }
}

void free_prio_queue(PrioQueue *q) {
  for(int i = 0; i != q->m_num_levels; i++)
    free(q->m_levels[i]);
  free(q->m_levels);
  free(q);
}
