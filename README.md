## 测试模板

```shell
# test
	- JUST_RUN_TEST
	- TEST
	- SKIP_TEST

	- ASSERT_***
	- EXPECT_***
# assert
	- 断言(不过并不会断掉程序)
# backtrace
	- 用于查看段错误(已经内置在main函数中)
# dbg
	- 用于测试简单的程序
# log
	- LOG("HELLO", "WORLD: ", 1) << std::endl;
	- LOG_FILE(ofs, "HELLO", "WORLD: ", 1) << std::endl;
```

### !注意

```txt
1. 使用多个JUST_RUN_TEST时会编译不通过, 会报: multiple definition of `g_justtest_Test'
```
