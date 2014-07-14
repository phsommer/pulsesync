#ifndef LOGICAL_CLOCK_H
#define LOGICAL_CLOCK_H

typedef struct TableItem
{
  uint32_t localTime;
  int32_t offset; // globalTime - localTime
} TableItem;


#endif /* LOGICAL_CLOCK_H */
