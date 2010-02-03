/**
**************************************************************************************
*Founder Effect                                                                      *
*Copyright (C) 2006-2007, 5d Studios (www.5d-Studios.com)                            *
*                                                                                    *
*This library is free software; you can redistribute it and/or                       *
*modify it under the terms of the GNU Lesser General Public                          *
*License as published by the Free Software Foundation; either                        *
*version 2.1 of the License, or (at your option) any later version.                  *
*                                                                                    *
*This library is distributed in the hope that it will be useful,                     *
*but WITHOUT ANY WARRANTY; without even the implied warranty of                      *
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU                   *
*Lesser General Public License for more details.                                     *
*                                                                                    *
*You should have received a copy of the GNU Lesser General Public                    *
*License along with this library; if not, write to the Free Software                 *
*Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA      *
**************************************************************************************
*/

#ifndef STRINGUTIL_H_
#define STRINGUTIL_H_

#include <string>
#include <list>
#include <stdlib.h>

namespace util
{


typedef std::list< std::string > type_StringList;
/**
 * A class for easy parsing of strings
 * it tokenizes the string from a desired tokenizer
 */
class StringUtil
{
public:
    StringUtil(void);
    /** default delimeter is ' ' */
    StringUtil( std::string s);
    /** Set a custom delimeter */
    StringUtil( std::string s, char delim);

    /** Test to see if there is a next token */
    bool HasNext() { return !m_tokens.empty(); };
    /** Get the next token */
    std::string GetNext() {
        std::string next = m_tokens.front();
        m_tokens.pop_front();
        return next;
    };

    /** Amount of tokens */
    size_t Size() { return m_tokens.size(); };

    /*------------------------------------------
        Util functions - probably shouldn't be here
        but owell.
    -------------------------------------------*/

    template<class T> std::string ToString( T val );


    static float ToFloat(const std::string &str) {
        return atof( str.c_str() );
    };

    static int ToInt( const std::string &str) {
        return atoi( str.c_str() );
    };

    static float ToFloat(char *str) {
        return atof( str );
    };

    static int ToInt( char *str) {
        return atoi( str );
    };

    static bool IsNumber( std::string &str );

private:
    /** Parse out the tokens */
    type_StringList Parse( std::string &s, char delim );
    // the current string
    std::string m_str;
    // tokens in the string
    type_StringList m_tokens;
public:
    virtual ~StringUtil(void);
};

}

#endif /* STRINGUTIL_H_ */
