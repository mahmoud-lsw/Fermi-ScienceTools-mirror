/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * $Id: MacOSDefs.hpp 568078 2007-08-21 11:43:25Z amassari $
 */


#ifndef MACOS_DEFS_HPP
#define MACOS_DEFS_HPP

// ---------------------------------------------------------------------------
//  MacOS runs in big endian mode (PPC) or little endian (intel).
// ---------------------------------------------------------------------------
// This had to be changed to work on intel Macs but is not needed with
// xerces version 3 and up.  jcv

#if defined(__BIG_ENDIAN__) || defined(__LITTLE_ENDIAN__)
#if __BIG_ENDIAN__
#define ENDIANMODE_BIG
#else
#define ENDIANMODE_LITTLE
#endif
#else
#define ENDIANMODE_BIG
#endif

// ---------------------------------------------------------------------------
//  Define all the required platform types
//
//	FileHandle is a pointer to XMLMacAbstractFile. Due to namespace
//	declaration issues, it is declared here as a void*.
// ---------------------------------------------------------------------------
typedef void*   FileHandle;

#endif
