// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             //  从 Windows 头文件中排除极少使用的信息
// Windows 头文件:
#include <windows.h>

// C 运行时头文件
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>
#include <iostream>
#include <sstream>

#define ENABLE_PARALLEL	1 // 多线程并行开关

#define NAME_STRING_PLATFORM_1	"Advanced Micro Devices, Inc."
#define NAME_STRING_PLATFORM	"Intel(R) OpenCL"
#define NAME_STRING_PLATFORM_3	"AMD Accelerated Parallel Processing"
#define NAME_STRING_PLATFORM_4	"NVIDIA CUDA"

#define  LocalWorkX		8
#define  LocalWorkY		8
#define  LocalWorkSizeDef	1//局部线程块采用默认数目

// TODO: 在此处引用程序需要的其他头文件
