/* Copyright (c) 2009, Distributed Computing Group (DCG), ETH Zurich.
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
* 
*/

interface TimeSync<clock_precision>
{

  /**
    * Return the current local time of the node. 
    * @return 	local timestamp
    */
  command uint32_t getLocalTime();
  
  
  /** Return the current global time of the node. This time may only be valid if the node
   *  is synchronized with the rest of the network. In this case the isSynced() method 
   *  will return TRUE. Otherwise the value returned by this method may be the same as the
   *  local time or it may have an arbitrary offset.
   * @return	global timestamp
   */
  command uint32_t getGlobalTime();
  
  
  /**
   * Convert a given local timestamp to the corresponding global timestamp.
   * @param 'uint32_t *localtime'	local timestamp
   * @return 	global timestamp
   */
  command void convertToGlobalTime(uint32_t *localtime);
  
  
  /**
   * Convert a given global timestamp to the corresponding local timestamp.
   * @param 'uint32_t *globaltime'	global timestamp
   * @return	local timestamp
   */
  command void convertToLocalTime(uint32_t *globaltime);
  
  
  /**
   * Returns TRUE if the node is synchronized (within certain bounds) with its neighboring nodes.
   * @return	whether node is synchronized
   */
  command bool isSynced();
  command float getSkew();


}

