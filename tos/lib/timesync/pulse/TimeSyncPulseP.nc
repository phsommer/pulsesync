/* Copyright (c) 2010, Distributed Computing Group (DCG), ETH Zurich.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *  3. Neither the name of the copyright holders nor the names of
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
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
#include "AM.h"
#include <stdint.h>
#include <math.h>

#ifdef TERMINAL_ENABLED
#include "printf.h"
#endif



module TimeSyncPulseP {
  provides
  {
    interface Init;
    interface StdControl; 

    interface TimeSync<PULSESYNC_CLOCK_TYPE>;
    interface TimeSyncInfo;
  }
  
  uses
  {
    
    interface LocalTime<PULSESYNC_CLOCK_TYPE>;
    interface LogicalClock as LogicalClock;
    interface LogicalClock as HardwareClock;
    
    interface Timer<TMilli> as BeaconTimer;
    interface Timer<TMilli> as SendTimer;
    
    interface SplitControl as RadioControl;

    interface Leds;
    interface Random;
    
    interface TimeSyncPacket<TRadio, uint32_t>;
    interface TimeSyncAMSend<TRadio, uint32_t>;
    interface PacketTimeStamp<TMilli, uint32_t> as PacketTimeStampMilli;
	
    interface Packet;
    interface AMPacket;
    interface Receive;

       
  }
}

implementation
{

  bool m_running = FALSE;
  bool isForwarding = FALSE;
  bool m_radioOn = FALSE;

  am_addr_t m_root;
  am_addr_t m_parent;

  timesync_point_t m_sync_point;  
  float m_init_skew;
  bool m_synchronized;
  uint16_t m_seqno;
  uint8_t m_hops;
  uint8_t m_heartbeat;
  
  /** message buffers **/
  message_t sync_msg;
  message_t* msg_ptr = &sync_msg;

  task void sendMsg();


  /** helper methods for timesync **/
  int32_t getError(uint32_t offset) {
    if (offset <= INT32_MAX) return offset;
    else return -(UINT32_MAX - offset);
  }
  
  bool isRoot() {
    return (m_root==TOS_NODE_ID);
  }
  
  void setRoot(am_addr_t root, am_addr_t parent, uint8_t hops, uint16_t seq) {

    if (root!=m_root) {
      call LogicalClock.reset();
      call HardwareClock.reset();
    }

    m_root = root; 
    m_hops = hops; 
    m_parent = parent;
    m_seqno = seq;
    
#ifdef TERMINAL_ENABLED 
    printf("setroot: %u %u %u %u\n", root, parent, hops, seq);
    printfflush();
#endif

    // reset heartbeat counter
    m_heartbeat = 0;
  }

  void sendBeacon(uint32_t backoff) {

    isForwarding = TRUE;

    if (backoff==0) {
      // send now
      call SendTimer.stop();
      post sendMsg();
    } else {
      call SendTimer.startOneShot(backoff);
    }
  }  
  
  task void updateClock() {
    
    int32_t error;
    uint32_t global = m_sync_point.localtime;
    signal TimeSyncInfo.updateReceived(m_parent, m_seqno, m_sync_point.localtime, m_sync_point.globaltime);


#ifdef PULSE_DRIFT_CORRECTION
    if (m_seqno % 2 == 0) {

      // beacons with even sequence numbers -> estimate local hardware clock drift
      call HardwareClock.addSyncEntry(m_sync_point.localtime, m_sync_point.globaltime);
 
    } else {
#endif
      // beacons with odd sequence numbers -> estimate offset/drift of logical clock

      call LogicalClock.convertToGlobalTime(&global);
   
      // check for consistency
      error = getError(m_sync_point.globaltime - global);
   
      if (m_synchronized && error>PULSE_MAX_ERROR_THRESHOLD) {

        call LogicalClock.reset();
        call HardwareClock.reset();

        #ifdef TERMINAL_ENABLED       
        printf("reset clocks: %ld\n", error);
        printfflush();
        #endif
      }

     // update regression table
     call LogicalClock.addSyncEntry(m_sync_point.localtime, m_sync_point.globaltime);
     m_synchronized = call LogicalClock.isSynchronized();

#ifdef PULSE_DRIFT_CORRECTION
    }
#endif
  } 
   

  task void sendMsg() {
    
    PulseSyncMsg_t* syncMsg = call Packet.getPayload(msg_ptr, sizeof(PulseSyncMsg_t));

    uint32_t localTime;
    int32_t offset = 0;
    
    correction_t c;
    c.f = 0.0f;
    syncMsg->correction = c.i;

    syncMsg->root = m_root;
    syncMsg->seqno = m_seqno;
    syncMsg->hops = m_hops;
    syncMsg->parent = m_parent;


    localTime = call LocalTime.get();

    syncMsg->localTime = localTime;

    // forward synchronization beacon
    if (!isRoot()) {
 
      double i_part, f_part;
      double corr = m_sync_point.correction.f;

      offset = m_sync_point.globaltime - m_sync_point.localtime;

      // add sub-tick part of message delay
      #ifdef PULSE_MESSAGE_DELAY_CONSTANT
      corr += PULSE_MESSAGE_DELAY_CONSTANT;
      #endif


      #ifdef PULSE_DRIFT_CORRECTION
      if (m_seqno % 2) {

        // correct local hardware clock drift for beacons with odd sequence numbers

        corr += call HardwareClock.getSkew()*((localTime - m_sync_point.localtime) + (196+8)); // add TX/RX_SFD_DELAY

      }
      #endif

      // split correction in integer and float parts
      f_part = corr;     
      i_part = rint(corr);
      f_part -= i_part;

      // add integer part to offset
      offset += i_part;

      c.f = f_part;
      syncMsg->correction = c.i;

    } 
	
    syncMsg->globalTime = localTime + offset;


    if (call TimeSyncAMSend.send(AM_BROADCAST_ADDR, msg_ptr, sizeof(PulseSyncMsg_t), localTime)!=SUCCESS) {
      isForwarding = FALSE;

      printf("sendFailed: radio: %u\n", m_radioOn);
      printfflush();

    }	
    

  }

  command error_t Init.init() {
    return SUCCESS;
  }

   
  /************ Interface StdControl **************************/
  command error_t StdControl.start()
  {
    error_t error;
    m_running = TRUE; 

    setRoot(AM_BROADCAST_ADDR, AM_BROADCAST_ADDR, 0xFF, 0);

    m_synchronized = FALSE;
    
    // reset clocks
    call LogicalClock.reset();
    call HardwareClock.reset();
    
    // start periodic beacon timer
    call BeaconTimer.startPeriodic(1024*PULSE_SYNC_PERIOD);

#ifdef TERMINAL_ENABLED
    printf("PulseSync started: %u\n", TOS_NODE_ID);
    printfflush();
#endif
 
    error = call RadioControl.start();
    if (error==EALREADY) {
      m_radioOn = TRUE;
#ifdef TERMINAL_ENABLED
      printf("radio already started\n");
      printfflush();
#endif
    }


    return SUCCESS;
  }

  command error_t StdControl.stop()
  {
    // stop timers
    call SendTimer.stop();
    call BeaconTimer.stop();
    m_running = FALSE;

    call RadioControl.stop();
    return SUCCESS;
  }

   
  /************ Interface Receive **************************/
  event message_t* Receive.receive(message_t* msg, void* payload, uint8_t len) {
    
    am_addr_t source = call AMPacket.source(msg);
     

    if (!m_running) return msg;

    if (call TimeSyncPacket.isValid(msg)==TRUE) {

      PulseSyncMsg_t* pulse = (PulseSyncMsg_t*)payload;
      uint16_t backoff;

      if (pulse->root==m_root && (int16_t)(pulse->seqno - m_seqno)>0) {

        // set root, parent, hops, seqno
        setRoot(pulse->root, source, pulse->hops + 1, pulse->seqno);

        m_sync_point.localtime = call TimeSyncPacket.eventTime(msg);
        m_sync_point.globaltime = pulse->globalTime;
        m_sync_point.correction.i = pulse->correction;
        
 
        backoff = PULSE_WAIT_BEFORE_SEND_INTERVAL;
        sendBeacon(backoff);

#ifdef TERMINAL_ENABLED
        //printf("update: local=%lu, global=%lu\n", m_sync_point.localtime, m_sync_point.globaltime);
        //printfflush();
#endif
                 
      } else if (pulse->root<m_root && pulse->root<TOS_NODE_ID) {

        // new root node
	setRoot(pulse->root, source, pulse->hops + 1, pulse->seqno);

        m_sync_point.localtime = call TimeSyncPacket.eventTime(msg);
        m_sync_point.globaltime = pulse->globalTime;
	m_sync_point.correction.i = pulse->correction;	

#ifdef TERMINAL_ENABLED
        printf("new root: %u\n", m_root);
	printfflush();
#endif
        
        // send immediately and stay awake
        sendBeacon(0);

      } else if (pulse->root>m_root) {
        // advertisement for a root node with a higher id than our current root
        // send immediately
        sendBeacon(0);

#ifdef TERMINAL_ENABLED
        printf("bogus root: %u\n", pulse->root);
	printfflush();
#endif
                
      } else if (pulse->root>TOS_NODE_ID) {
	// we are a better candidate for the new root but still waiting
	setRoot(TOS_NODE_ID, TOS_NODE_ID, 0, pulse->seqno);
        sendBeacon(0);
      } else {
        // ignore beacon
        
      }


    } else {
#ifdef TERMINAL_ENABLED
      printf("timestamp invalid\n");
      printfflush();
#endif  
    }

    return msg;
  }
     
  /************ Interface timeSyncAMSend **************************/
  event void TimeSyncAMSend.sendDone(message_t* msg, error_t error) {

      msg_ptr = msg;
      isForwarding = FALSE;

      // start overhearing timer
      //call OverhearTimer.startOneShot(PULSE_AFTER_SLOT_INTERVAL);

      if (!isRoot()) {
        post updateClock();
      }

  }
  
  
 
  /************ Interface Timer **************************/
  event void BeaconTimer.fired() {

    if (isRoot()) {

      // increase sequence number
      m_seqno++;  

      if (m_radioOn) {
        // radio started, send message immediately
        sendBeacon(0);
      } else {
        // start radio
        error_t error = call RadioControl.start();
	if (error!=SUCCESS) { 
#ifdef TERMINAL_ENABLED
          printf("Failed to start radio: %u\n", error);
	  printfflush();
#endif
        }
      }

    } else if (++m_heartbeat==PULSE_IGNORE_ROOT_ROUNDS) {
      // declare myself the new root node
      setRoot(TOS_NODE_ID, TOS_NODE_ID, 0, m_seqno);
 
#ifdef TERMINAL_ENABLED
      printf("declared root\n");
      printfflush();
#endif
      
      // radio should be enabled here since we are listening for other beacons
      sendBeacon(0);
    }
  }
  
  event void SendTimer.fired() {
     post sendMsg();
  }


  /************ Interface RadioControl ***********************/
  event void RadioControl.startDone(error_t error){

    if (error==SUCCESS) {
        // radio started
        m_radioOn = TRUE;
      	if (m_running && isRoot()) {
          post sendMsg();
        }
    } else {
#ifdef TERMINAL_ENABLED
      printf("Failed to start radio: %u\n", error);
      printfflush();
#endif
    }
  }

  event void RadioControl.stopDone(error_t error) {

    if (error==SUCCESS) {
        // radio stopped
        m_radioOn = FALSE;
    } else {
#ifdef TERMINAL_ENABLED
      printf("Failed to stop radio: %u\n", error);
      printfflush();
#endif
    }
  }

 
  /************ Interface TimeSync **************************/
  command uint32_t TimeSync.getLocalTime()
  {
    return call LocalTime.get();
  }

  command uint32_t TimeSync.getGlobalTime()
  {
    uint32_t time = call TimeSync.getLocalTime();
    call LogicalClock.convertToGlobalTime(&time);
    return time;
  }

  command void TimeSync.convertToGlobalTime(uint32_t *time)
  {
    call LogicalClock.convertToGlobalTime(time);
  }

  command void TimeSync.convertToLocalTime(uint32_t *time)
  {
    call LogicalClock.convertToLocalTime(time);
  }

  command bool TimeSync.isSynced(){
    return m_synchronized || (m_root == TOS_NODE_ID);  
  }
  
  command float TimeSync.getSkew() {
    return call LogicalClock.getSkew();
  }

  command am_addr_t TimeSyncInfo.getRoot() {
    return m_root; 
  }  

  command am_addr_t TimeSyncInfo.getParent() {
    return m_parent; 
  }  

  command bool TimeSyncInfo.isStarted() {
    return m_running; 
  }  
  
}
