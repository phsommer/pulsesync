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

configuration TimeSyncPulseC {
  provides interface Init;
  provides interface StdControl;
  provides interface TimeSync<PULSESYNC_CLOCK_TYPE>;
  provides interface TimeSyncInfo;
}


implementation {
 
  components TimeSyncPulseP as TimeSyncP;
  Init = TimeSyncP;
  StdControl = TimeSyncP;
  TimeSync = TimeSyncP;
  TimeSyncInfo = TimeSyncP;
  
  components new TimerMilliC() as BeaconTimer;
  TimeSyncP.BeaconTimer -> BeaconTimer;
  
  components new TimerMilliC() as SendTimer;
  TimeSyncP.SendTimer -> SendTimer;

  components ActiveMessageC, TimeSyncMessageC;
  TimeSyncP.Packet -> ActiveMessageC;
  TimeSyncP.AMPacket -> ActiveMessageC;
  TimeSyncP.Receive ->  TimeSyncMessageC.Receive[TIMESYNC_AM_PULSE];
  TimeSyncP.RadioControl -> ActiveMessageC;
 
#ifdef CLOCK_PRECISION_MILLI
  TimeSyncP.TimeSyncAMSend -> TimeSyncMessageC.TimeSyncAMSendMilli[TIMESYNC_AM_PULSE];
  TimeSyncP.TimeSyncPacket -> TimeSyncMessageC.TimeSyncPacketMilli;
  components LocalTimeMilliC;
  TimeSyncP.LocalTime -> LocalTimeMilliC;
#endif

#ifdef CLOCK_PRECISION_RADIO
  TimeSyncP.TimeSyncAMSend -> TimeSyncMessageC.TimeSyncAMSendRadio[TIMESYNC_AM_PULSE];
  TimeSyncP.TimeSyncPacket -> TimeSyncMessageC.TimeSyncPacketRadio;
  TimeSyncP.PacketTimeStampMilli -> TimeSyncMessageC.PacketTimeStampMilli;

  #if defined (PLATFORM_RCB128RFA1)
  components LocalTime62khzC as LocalTime;
  #elif defined (PLATFORM_TELOSB)
  components LocalTime32khzC as LocalTime;
  #elif defined (PLATFORM_IRIS) || defined(PLATFORM_OPAL)
  components LocalTimeMicroC as LocalTime;
  #else
  #error "Unsupported platform"
  #endif
  TimeSyncP.LocalTime -> LocalTime;
#endif

  components new LogicalClockC(PULSE_CLOCK_TABLE_SIZE) as LogicalClock, new LogicalClockC(PULSE_CLOCK_TABLE_SIZE) as HardwareClock;
  TimeSyncP.LogicalClock -> LogicalClock;  
  TimeSyncP.HardwareClock -> HardwareClock; 
    
  components LedsC;
  TimeSyncP.Leds -> LedsC;
  
  components RandomC;
  TimeSyncP.Random -> RandomC;

	
}
