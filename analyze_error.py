#!/usr/bin/env python3
"""
编译错误自动归因与修复建议脚本
功能：分析编译错误日志，提供针对性修复建议
"""

import re
import sys

# 错误-解决方案映射表
error_solutions = {
    # 缺失MPI头文件
    r"fatal error: mpi.h: No such file or directory": {
        "原因": "缺少MPI开发库",
        "解决方案": [
            "Ubuntu/Debian: sudo apt-get install libopenmpi-dev",
            "CentOS/RHEL: sudo yum install openmpi-devel",
            "macOS: brew install open-mpi",
            "编译时确保启用MPI: cmake -DMPI_FOUND=ON .."
        ]
    },
    # 缺失Boost库
    r"fatal error: boost/.*: No such file or directory": {
        "原因": "缺少Boost开发库",
        "解决方案": [
            "Ubuntu/Debian: sudo apt-get install libboost-all-dev",
            "CentOS/RHEL: sudo yum install boost-devel",
            "macOS: brew install boost"
        ]
    },
    # 缺失Eigen库
    r"fatal error: Eigen/.*: No such file or directory": {
        "原因": "缺少Eigen开发库",
        "解决方案": [
            "Ubuntu/Debian: sudo apt-get install libeigen3-dev",
            "CentOS/RHEL: sudo yum install eigen3-devel",
            "macOS: brew install eigen"
        ]
    },
    # 缺失HDF5库
    r"fatal error: hdf5/.*: No such file or directory": {
        "原因": "缺少HDF5开发库",
        "解决方案": [
            "Ubuntu/Debian: sudo apt-get install libhdf5-dev",
            "CentOS/RHEL: sudo yum install hdf5-devel",
            "macOS: brew install hdf5"
        ]
    },
    # 语法错误
    r"error: .*syntax error": {
        "原因": "C++语法错误",
        "解决方案": [
            "检查错误行附近的语法",
            "确保使用正确的C++标准（项目要求C++14）"
        ]
    },
    # 未定义引用（MPI）
    r"undefined reference to.*MPI_": {
        "原因": "MPI库链接错误",
        "解决方案": [
            "确保编译时启用了MPI支持",
            "检查MPI库路径是否正确",
            "尝试重新运行cmake配置"
        ]
    },
    # 未定义引用（Boost）
    r"undefined reference to.*boost::": {
        "原因": "Boost库链接错误",
        "解决方案": [
            "确保Boost库已正确安装",
            "检查Boost库版本是否兼容",
            "尝试重新运行cmake配置"
        ]
    },
    # CMake配置错误
    r"CMake Error at.*": {
        "原因": "CMake配置错误",
        "解决方案": [
            "删除build目录并重新配置: rm -rf build && mkdir build && cd build && cmake ..",
            "检查CMake版本是否符合要求（项目要求>=3.12）",
            "检查是否缺少必要的依赖包"
        ]
    },
    # 权限错误
    r"permission denied": {
        "原因": "权限不足",
        "解决方案": [
            "确保有足够的权限执行编译操作",
            "尝试使用sudo（谨慎使用）"
        ]
    },
    # 内存不足
    r"collect2: error: ld returned 1 exit status.*out of memory": {
        "原因": "编译时内存不足",
        "解决方案": [
            "减少并行编译的线程数: make -j$(nproc-1)",
            "增加系统内存或使用交换空间"
        ]
    }
}

def analyze_error_log(log_content):
    """分析错误日志并返回修复建议"""
    suggestions = []
    
    for error_pattern, solution_info in error_solutions.items():
        if re.search(error_pattern, log_content):
            suggestions.append({
                "原因": solution_info["原因"],
                "解决方案": solution_info["解决方案"]
            })
    
    return suggestions

def main():
    if len(sys.argv) != 2:
        print("用法: python3 analyze_error.py <错误日志文件>")
        print("或: 编译命令 2>&1 | python3 analyze_error.py -")
        sys.exit(1)
    
    log_file = sys.argv[1]
    
    if log_file == "-":
        # 从标准输入读取
        log_content = sys.stdin.read()
    else:
        # 从文件读取
        try:
            with open(log_file, 'r') as f:
                log_content = f.read()
        except FileNotFoundError:
            print(f"错误: 文件 {log_file} 不存在")
            sys.exit(1)
    
    suggestions = analyze_error_log(log_content)
    
    if suggestions:
        print("=== 编译错误分析结果 ===")
        for i, suggestion in enumerate(suggestions, 1):
            print(f"\n问题 {i}: {suggestion['原因']}")
            print("解决方案:")
            for sol in suggestion['解决方案']:
                print(f"  - {sol}")
    else:
        print("未识别到已知的编译错误模式。")
        print("建议检查:")
        print("1. 依赖库是否完整安装")
        print("2. CMake配置是否正确")
        print("3. 代码是否符合C++14标准")

if __name__ == "__main__":
    main()
