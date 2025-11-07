#!/bin/bash

# 编译环境快照与回滚脚本 (Linux/macOS)
# 功能：备份和恢复编译环境的关键配置

ACTION="status"
SNAPSHOT_FILE="env_snapshot.txt"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -a|--action)
            ACTION="$2"
            shift # past argument
            shift # past value
            ;;
        -f|--file)
            SNAPSHOT_FILE="$2"
            shift # past argument
            shift # past value
            ;;
        -h|--help)
            echo "用法: $0 [选项]"
            echo "选项:"
            echo "  -a, --action ACTION   操作类型: backup, restore, status (默认: status)"
            echo "  -f, --file FILE       快照文件名 (默认: env_snapshot.txt)"
            echo "  -h, --help            显示帮助信息"
            exit 0
            ;;
        *)
            echo "未知选项: $1"
            echo "使用 -h 查看帮助"
            exit 1
            ;;
    esac
done

backup_environment() {
    # 备份环境配置
    echo "正在创建编译环境快照..."
    
    cat > "$SNAPSHOT_FILE" <<EOF
=== 系统信息 ===
操作系统: $(uname -a)
当前用户: $(whoami)
当前目录: $(pwd)
日期: $(date)

=== 环境变量 ===
CXX_COMPILER: ${CXX_COMPILER:-未设置}
CC_COMPILER: ${CC_COMPILER:-未设置}
CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH:-未设置}
BOOST_ROOT: ${BOOST_ROOT:-未设置}
EIGEN3_ROOT: ${EIGEN3_ROOT:-未设置}
MPI_ROOT: ${MPI_ROOT:-未设置}
HDF5_ROOT: ${HDF5_ROOT:-未设置}
PATH: $PATH
LD_LIBRARY_PATH: ${LD_LIBRARY_PATH:-未设置}

=== MPI相关库 ===
$(ls /usr/lib/*mpi* 2>/dev/null || echo "未找到MPI库")
$(ls /usr/include/*mpi* 2>/dev/null || echo "未找到MPI头文件")

=== CMake信息 ===
$(cmake --version 2>/dev/null || echo "CMake未安装")

=== 编译器信息 ===
CXX编译器: ${CXX_COMPILER:-默认编译器}
$(if [ -n "$CXX_COMPILER" ]; then $CXX_COMPILER --version 2>&1 | head -1; else g++ --version 2>&1 | head -1; fi)
CC编译器: ${CC_COMPILER:-默认编译器}
$(if [ -n "$CC_COMPILER" ]; then $CC_COMPILER --version 2>&1 | head -1; else gcc --version 2>&1 | head -1; fi)
EOF
    
    echo "环境快照已保存到: $SNAPSHOT_FILE"
}

restore_environment() {
    # 恢复环境配置
    if [ ! -f "$SNAPSHOT_FILE" ]; then
        echo "错误: 快照文件不存在: $SNAPSHOT_FILE"
        return 1
    fi
    
    echo "正在恢复编译环境..."
    
    # 恢复环境变量
    while IFS=: read -r var_name var_value; do
        var_name=$(echo "$var_name" | xargs)
        var_value=$(echo "${var_value:-}" | xargs)
        
        case "$var_name" in
            CXX_COMPILER|CC_COMPILER|CMAKE_PREFIX_PATH|BOOST_ROOT|EIGEN3_ROOT|MPI_ROOT|HDF5_ROOT)
                if [ -n "$var_value" ] && [ "$var_value" != "未设置" ]; then
                    export "$var_name=$var_value"
                    echo "已恢复环境变量: $var_name=$var_value"
                fi
                ;;
        esac
    done < <(grep -E "^(CXX_COMPILER|CC_COMPILER|CMAKE_PREFIX_PATH|BOOST_ROOT|EIGEN3_ROOT|MPI_ROOT|HDF5_ROOT):" "$SNAPSHOT_FILE")
    
    echo "环境恢复完成！"
    echo "注意: 已安装的依赖包无法自动回滚，请手动卸载新增的包"
}

show_status() {
    # 显示当前环境状态
    echo "=== 当前编译环境状态 ==="
    
    # 显示环境变量
    env_vars=("CXX_COMPILER" "CC_COMPILER" "CMAKE_PREFIX_PATH" "BOOST_ROOT" "EIGEN3_ROOT" "MPI_ROOT" "HDF5_ROOT")
    for var in "${env_vars[@]}"; do
        value=${!var:-未设置}
        echo "$var: $value"
    done
    
    # 显示CMake版本
    if command -v cmake &> /dev/null; then
        echo "CMake版本: $(cmake --version | head -1)"
    else
        echo "CMake: 未安装"
    fi
    
    # 显示快照文件状态
    if [ -f "$SNAPSHOT_FILE" ]; then
        snapshot_time=$(stat -c "%y" "$SNAPSHOT_FILE" 2>/dev/null || stat -f "%Sm" "$SNAPSHOT_FILE")
        echo "最新快照: $SNAPSHOT_FILE (创建于: $snapshot_time)"
    else
        echo "快照文件: 不存在"
    fi
}

# 主逻辑
case "$ACTION" in
    backup)
        backup_environment
        ;;
    restore)
        restore_environment
        ;;
    status)
        show_status
        ;;
    *)
        echo "无效的操作: $ACTION"
        echo "可用操作: backup, restore, status"
        exit 1
        ;;
esac
