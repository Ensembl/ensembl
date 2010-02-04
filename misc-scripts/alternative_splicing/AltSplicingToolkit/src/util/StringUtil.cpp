#include <sstream>

#include "StringUtil.h"

namespace util
{


StringUtil::StringUtil(void)
{
}

/** default delimeter is ' ' */
StringUtil::StringUtil(std::string s, char delim)
{
    m_tokens = Parse( s, delim );
}

/** Set a custom delimeter */
StringUtil::StringUtil(std::string s)
{
    m_tokens = Parse( s, ' ');
}

StringUtil::~StringUtil(void)
{
}

/** Parse out the tokens */
type_StringList StringUtil::Parse( std::string &str,  char delim )
{
    type_StringList tokens;
    // initial spot for substr
    unsigned int j = 0, i = 0;

    int size = 0;
    // search through the whole string to
    // find tokens
    while(i < str.length() && j < str.length() )
    {
        size = 0;
        // eat white space
        while( str[j] == delim && j < str.length() ) {
            j++;
            if ( j < str.length() )
                i = j;
        }

        // if we find a space
        while( str[i] != delim && i < str.length() )
        { 
            i++;
            size++;
        }
        if ( j <= i-1 )
            // add this token to the stack
            tokens.push_back( str.substr(j, size) );
        if ( j >= str.length() )
            return tokens;
        j = i;
    }
    return tokens;
}


template<class T>
std::string StringUtil::ToString( T val )
{
    std::ostringstream sstream;
    sstream << val;
    return sstream.str();
}


// helper function
bool IsNumber_c( char c )
{
    // a number
    char nums[] = { '0', '1', '2', '3', '4',
                    '5', '6', '7', '8', '9', '.'
    };

    // search the array
    for ( int i = 0; i < 11; i++ )
    {
        // if we have a match
        // it is a number
        if ( c == nums[i] )
            return true;
    }
    return false;
}

/** Check for a number 
    Note very fragile - TODO add a better check for decimals
 */
bool StringUtil::IsNumber( std::string &str )
{
    for ( unsigned int i = 0; i < str.length(); i++ )
    {
        if ( !IsNumber_c( str[i] ) )
            return false;
    }
    return true;
}

}
