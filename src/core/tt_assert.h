#ifndef __ASSERT_H__
#define __ASSERT_H__

#include<iostream>
#define tt_assert(expr, msg) \
    if(!(expr))  std::cout << "Assertion failed: [" << __FILE__ << " line " << __LINE__ << "] : " << #expr <<  " ( " << msg << " )" << std::endl;
       

#endif // !__ASSERT_H__