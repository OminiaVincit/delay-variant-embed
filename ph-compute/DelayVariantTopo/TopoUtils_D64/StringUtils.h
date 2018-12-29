#pragma once
#include "TopoUtils_D64.h"
#include <string>
#include <comutil.h>
#include <strsafe.h>
#include <sstream>

namespace NStringUtil {
#ifndef UNICODE
    typedef std::string tstring;
    typedef char tchar_t;
#else
    typedef std::wstring tstring;
    typedef wchar_t tchar_t;
#endif

    TopoUtils_D64_API std::wstring _s2w(const std::string &src, UINT codepage = ::GetACP());
    TopoUtils_D64_API std::string _w2s(const std::wstring &src, UINT codepage = ::GetACP());
};