#include"core/tt_test.h"
#include"core/tt_backtrace.h"

int main()
{
    REGISTER_SEGFAULT_HANDLER
    
    return RunAllTests();
}