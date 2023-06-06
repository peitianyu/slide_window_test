#ifndef __ASSERT_H__
#define __ASSERT_H__

#include<iostream>
#define tt_assert(expr) \
    if(!(expr))  std::cout << "Assertion failed: [" << __FILE__ << " line " << __LINE__ << "] : " << #expr <<  std::endl; 
       

#endif // !__ASSERT_H__