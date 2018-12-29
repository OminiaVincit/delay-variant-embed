#pragma once

// need to use the latest VC++ runtime libraries
#if !defined(_BIND_TO_CURRENT_VCLIBS_VERSION) || _BIND_TO_CURRENT_VCLIBS_VERSION != 1
#error "_BIND_TO_CURRENT_VCLIBS_VERSION must be defined as 1. Use Common.vsprops"
#endif

// for Win32 support WinXP SP2 or later
// for x64 support Vista/Server2008 or later
// See http://msdn.microsoft.com/en-us/library/aa383745(VS.85).aspx for more information
// 
// If you need to override WINVER       define __OVERRIDE_WINVER__      before including this header (typically inside stdafx.h)
// If you need to override _WIN32_WINNT define __OVERRIDE_WIN32_WINNT__ before including this header (typically inside stdafx.h)
// If you need to override _WIN32_IE    define __OVERRIDE_WIN32_IE__    before including this header (typically inside stdafx.h)
// 
#ifndef WINVER
#define WINVER _WIN32_WINNT_WIN7
#else
#error "WINVER must not be defined elsewhere ! Please use the one defined here."
#endif

#ifndef _WIN32_WINNT
#define _WIN32_WINNT _WIN32_WINNT_WIN7
#else
#error "_WIN32_WINNT must not be defined elsewhere ! Please use the one defined here."
#endif

#ifndef _WIN32_IE
#define _WIN32_IE _WIN32_IE_IE80
#else
#error "_WIN32_IE must not be defined elsewhere ! Please use the one defined here."
#endif

#define _USE_MATH_DEFINES
#define STRSAFE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#define _SCL_SECURE_NO_DEPRECATE
#define _ATL_SECURE_NO_DEPRECATE
#define _AFX_SECURE_NO_DEPRECATE
#if !defined(_SECURE_SCL) || _SECURE_SCL != 0
#error "_SECURE_SCL must be defined as 0. Use Common.vsprops"
#endif	// #if not defined(_SECURE_SCL) or _SECURE_SCL != 0
#ifndef _WIN64
#define _USE_32BIT_TIME_T
#endif

// Ignoring warning on warning C4996: 'MBCS_Support_Deprecated_In_MFC': MBCS support in MFC is deprecated and may be removed in a future version of MFC
// See http://stackoverflow.com/questions/24244544/visual-studio-2013-doesnt-ignore-disabled-warnings
#define NO_WARN_MBCS_MFC_DEPRECATION

#include <list>
#include <vector>
#include <deque>
#include <map>
// HashMap is no longer supported in VS2015
// #include <hash_map>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>

#pragma warning(push)
#pragma warning(disable : 4244 4267)
#include <boost/format.hpp>
#include <boost/utility.hpp>
#include <boost/cstdint.hpp>
#pragma warning(pop)

#include <boost/cstdint.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/optional.hpp>
#include <boost/version.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/preprocessor.hpp>

#if BOOST_VERSION >= 103200

	#include <boost/bind.hpp>
	#include <boost/format.hpp>

#endif


