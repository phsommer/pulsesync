/* Copyright (c) 2010, Distributed Computing Group (DCG), ETH Zurich.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*  1. Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*  2. Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in the
*     documentation and/or other materials provided with the distribution.
*  3. Neither the name of the copyright holders nor the names of
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS `AS IS'
*  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
*  ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
*  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, LOSS OF USE, DATA,
*  OR PROFITS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
*  THE POSSIBILITY OF SUCH DAMAGE.
*
*  @author Philipp Sommer <sommer@tik.ee.ethz.ch>
* 
*/

#include "TimeSyncPulse.h"
#include "LogicalClock.h"


generic module LogicalClockC(uint8_t tableSize) {
  provides interface LogicalClock;
}

implementation {


  /*  
      Based on the original FTSP implementation:

      We do linear regression from localTime to timeOffset (globalTime - localTime).
      This way we can keep the slope close to zero (ideally) and represent it
      as a float with high precision.

        timeOffset - offsetAverage = skew * (localTime - localAverage)
        timeOffset = offsetAverage + skew * (localTime - localAverage)
        globalTime = localTime + offsetAverage + skew * (localTime - localAverage)
  */
  
  TableItem table[tableSize];
  
  uint8_t next;
  uint8_t entries;

  float skew;
  uint32_t localAverage;
  int32_t offsetAverage;
 

  command void LogicalClock.reset() {
    next = 0;
    entries = 0;
    skew = 0.0f;
    localAverage = 0;
    offsetAverage = 0;
  }
  
  command void LogicalClock.addSyncEntry(uint32_t localTime, uint32_t globalTime) {
        
    int32_t localAverageRest = 0;
    int32_t offsetAverageRest = 0;
    int64_t localSum = 0;
    int64_t offsetSum = 0;

    uint8_t i;
        
    // insert entry into the circular buffer
    uint32_t newLocalAverage = table[next].localTime = localTime;
    int32_t newOffsetAverage = table[next].offset = globalTime - localTime;
    
    // update table size
    if (entries<tableSize) entries++;
    
#ifdef TERMINAL_ENABLED
    //printf("add: local=%lu, global=%lu\n", localTime, globalTime);
#endif

   /*
    We use a rough approximation first to avoid time overflow errors. The idea
    is that all times in the table should be relatively close to each other.
    */

    for (i=0; i<entries; i++) {
         
         
#ifdef TERMINAL_ENABLED
    //printf("Table: %u, local: %lu, offset: %ld\n", i, table[i].localTime, table[i].offset);       
#endif  
         
      /*
       This only works because C ISO 1999 defines the sign for modulo the same as for the Dividend!
      */ 
            	
       localSum += (int32_t)(table[i].localTime - newLocalAverage) / entries;
       localAverageRest += (table[i].localTime - newLocalAverage) % entries;
                
       offsetSum += (int32_t)(table[i].offset - newOffsetAverage) / entries;
       offsetAverageRest += (table[i].offset - newOffsetAverage) % entries;
     }

     newLocalAverage += localSum + ((localAverageRest + entries/2) / entries);
     newOffsetAverage += offsetSum + ((offsetAverageRest + entries/2) / entries);
        

     localSum = offsetSum = 0;
     for (i = 0; i < entries; ++i) {
       int32_t a = table[i].localTime - newLocalAverage;
       int32_t b = table[i].offset - newOffsetAverage;

       localSum += (int64_t)a * a;
       offsetSum += (int64_t)a * b;
     }

     if (localSum != 0) skew = (float)offsetSum / (float)localSum;
     else skew = 0.0f;

     offsetAverage = newOffsetAverage;
     localAverage = newLocalAverage;

#ifdef TERMINAL_ENABLED
     //printf("localAverage: %lu, offsetAverage: %ld\n", localAverage, offsetAverage);
     //printf("skew: %f\n", skew);
#endif

    // move pointer to the next free entry
    next = (next + 1) % tableSize;
    
  }
  
  command uint8_t LogicalClock.getEntries() {
    return entries;
  }
  
  command bool LogicalClock.isSynchronized() {
    return entries==tableSize;
  }
  
  command float LogicalClock.getSkew() {
    return skew;
  }
  

  command int32_t LogicalClock.getOffset() {
    return offsetAverage;
  }

  command void LogicalClock.convertToGlobalTime(uint32_t *time) {
    *time += offsetAverage + (int32_t)(skew * (int32_t)(*time - localAverage));
  }
  
  command void LogicalClock.convertToLocalTime(uint32_t *time) {
    uint32_t approxLocalTime = *time - offsetAverage;
    *time = approxLocalTime - (int32_t)(skew * (int32_t)(approxLocalTime - localAverage));
  }


}
