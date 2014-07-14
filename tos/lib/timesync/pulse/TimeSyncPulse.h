#ifndef TIMESYNC_PULSE_H
#define TIMESYNC_PULSE_H


#include "AM.h"

#define PULSE_SYNC_PERIOD 10L
#define PULSE_NEIGHBOR_TABLE_SIZE 16
#define PULSE_MAX_NEIGHBOR_AGE 3
#define PULSE_MAX_ERROR_THRESHOLD 1000
#define PULSE_IGNORE_ROOT_ROUNDS 5 // how many rounds to wait before declaring root


#ifndef PULSE_WAIT_BEFORE_SEND_INTERVAL
#define PULSE_WAIT_BEFORE_SEND_INTERVAL 0L
#endif

#define PULSE_PRE_SLOT_INTERVAL 10L
#define PULSE_AFTER_SLOT_INTERVAL 100L
#define PULSE_AWAKE_INTERVAL 3

#ifndef PULSE_CLOCK_TABLE_SIZE
#define PULSE_CLOCK_TABLE_SIZE 4
#endif

#define TIMESYNC_AM_PULSE 0x40


#ifdef CLOCK_PRECISION_RADIO
#define PULSESYNC_CLOCK_TYPE TRadio
#endif

#ifdef CLOCK_PRECISION_MILLI
#define PULSESYNC_CLOCK_TYPE TMilli
#endif

typedef union
{
  int32_t i;
  float f;
} correction_t;


/** Structure of the PulseSync Message **/
typedef nx_struct PulseSyncMsg {
  nx_uint8_t hops;
  nx_am_addr_t root;
  nx_uint16_t seqno;
  nx_am_addr_t parent;
  nx_uint32_t localTime;
  nx_uint32_t globalTime;
  nx_int32_t correction;
} PulseSyncMsg_t;



typedef struct timesync_point
{
  uint32_t localtime;
  uint32_t globaltime;
  correction_t correction;
} timesync_point_t;




//#define PULSESYNC_DEBUG

#endif /* TIME_SYNC_PULSE_H */
