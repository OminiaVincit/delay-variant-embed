#pragma once
#include "stdafx.h"
#include "StringUtils.h"

namespace NStringUtil {

    inline std::wstring stringToWString(const std::string &src, UINT codepage = ::GetACP()) {
        std::wstring ret;
        if (src.empty()) return ret;
        if (codepage == 0) {
            // convert 8bit -> 16bit
            ret.resize(src.size());
            std::string::const_iterator ps = src.begin();
            std::wstring::iterator pd = ret.begin();
            while (ps != src.end())
                *pd++ = (unsigned char)*ps++;
        }
        else {
            ret.resize(src.size());
            int sz = ::MultiByteToWideChar(codepage, 0, &src[0], (int)src.size(), &ret[0], (int)ret.size());
            if (sz == 0 && ::GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
                sz = ::MultiByteToWideChar(codepage, 0, &src[0], (int)src.size(), NULL, 0);
                ret.resize(sz);
                sz = ::MultiByteToWideChar(codepage, 0, &src[0], (int)src.size(), &ret[0], (int)ret.size());
            }
            ret.resize(sz);
        }
        return ret;
    };
    TopoUtils_D64_API std::wstring _s2w(const std::string &src, UINT codepage) { 
        return stringToWString(src, codepage); 
    };

    inline std::string  wstringToString(const std::wstring &src, UINT codepage = ::GetACP()) {
        std::string ret;
        if (src.empty()) return ret;
        if (codepage == 0) {
            // convert 16bit -> 8bit
            ret.resize(src.size());
            std::wstring::const_iterator ps = src.begin();
            std::string::iterator pd = ret.begin();
            while (ps != src.end())
                *pd++ = (unsigned char)*ps++;
        }
        else {
            ret.resize(src.size() * 2 + 256);
            int sz = ::WideCharToMultiByte(codepage, 0, &src[0], (int)src.size(), &ret[0], (int)ret.size(), NULL, NULL);
            if (sz == 0 && ::GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
                sz = ::WideCharToMultiByte(codepage, 0, &src[0], (int)src.size(), NULL, 0, NULL, NULL);
                ret.resize(sz);
                sz = ::WideCharToMultiByte(codepage, 0, &src[0], (int)src.size(), &ret[0], (int)ret.size(), NULL, NULL);
            }
            ret.resize(sz);
        }
        return ret;
    }
    TopoUtils_D64_API std::string _w2s(const std::wstring &src, UINT codepage) { 
        return wstringToString(src, codepage); 
    };

#ifndef UNICODE
    inline std::wstring tstringToWString(const tstring &src, UINT codepage = ::GetACP()) {
        return stringToWString(src, codepage);
    }

    inline tstring wstringToTString(const std::wstring &src, UINT codepage = ::GetACP()) {
        return wstringToString(src, codepage);
    }

    inline std::string tstringToString(const tstring &src, UINT codepage = ::GetACP()) {
        return src;
    }
    inline std::string tstringToString(tstring &&src, UINT codepage = ::GetACP()) {
        return std::move(src);
    }

    inline tstring stringToTString(const std::string &src, UINT codepage = ::GetACP()) {
        return src;
    }
    inline tstring stringToTString(std::string &&src, UINT codepage = ::GetACP()) {
        return std::move(src);
    }

    inline std::wstring _t2w(const tstring &src, UINT codepage = ::GetACP()) { return stringToWString(src, codepage); };
    inline tstring      _w2t(const std::wstring &src, UINT codepage = ::GetACP()) { return wstringToString(src, codepage); };
    inline std::string  _t2s(const tstring &src, UINT codepage = ::GetACP()) { return src; };
    inline std::string  _t2s(tstring &&src, UINT codepage = ::GetACP()) { return std::move(src); };
    inline tstring      _s2t(const std::string &src, UINT codepage = ::GetACP()) { return src; };
    inline tstring      _s2t(std::string &&src, UINT codepage = ::GetACP()) { return std::move(src); };
#else
    inline std::wstring tstringToWString(const tstring &src, UINT codepage = ::GetACP()) {
        return src;
    }
    inline std::wstring tstringToWString(tstring &&src, UINT codepage = ::GetACP()) {
        return std::move(src);
    }
    inline tstring wstringToTString(const std::wstring &src, UINT codepage = ::GetACP()) {
        return src;
    }
    inline tstring wstringToTString(std::wstring &&src, UINT codepage = ::GetACP()) {
        return std::move(src);
    }

    inline std::string tstringToString(const tstring &src, UINT codepage = ::GetACP()) {
        return wstringToString(src, codepage);
    }

    inline tstring stringToTString(const std::string &src, UINT codepage = ::GetACP()) {
        return stringToWString(src, codepage);
    }
    inline std::wstring _t2w(const tstring &src, UINT codepage = ::GetACP()) { return src; };
    inline std::wstring _t2w(tstring &&src, UINT codepage = ::GetACP()) { return std::move(src); };
    inline tstring      _w2t(const std::wstring &src, UINT codepage = ::GetACP()) { return src; };
    inline tstring      _w2t(std::wstring &&src, UINT codepage = ::GetACP()) { return std::move(src); };
    inline std::string  _t2s(const tstring &src, UINT codepage = ::GetACP()) { return wstringToString(src, codepage); };
    inline tstring      _s2t(const std::string &src, UINT codepage = ::GetACP()) { return stringToWString(src, codepage); };
#endif

};