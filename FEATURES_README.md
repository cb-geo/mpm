# 编译错误自动归因与修复建议

## 功能介绍
该功能可以自动分析编译或测试过程中产生的错误日志，识别常见的错误类型（如缺失依赖、语法错误、配置冲突等），并提供针对性的修复建议。

## 支持的错误类型
- 缺失MPI开发库
- 缺失Boost开发库  
- 缺失Eigen开发库
- 缺失HDF5开发库
- C++语法错误
- MPI库链接错误
- Boost库链接错误
- CMake配置错误
- 权限错误
- 内存不足错误

## 使用方法

### 1. 直接分析错误日志文件
```bash
python3 analyze_error.py <错误日志文件>
```

### 2. 实时分析编译输出
```bash
cmake --build build 2>&1 | python3 analyze_error.py -
```
或
```bash
make -j4 2>&1 | python3 analyze_error.py -
```

### 3. 示例输出
```
=== 编译错误分析结果 ===

问题 1: 缺少MPI开发库
解决方案:
  - Ubuntu/Debian: sudo apt-get install libopenmpi-dev
  - CentOS/RHEL: sudo yum install openmpi-devel
  - macOS: brew install open-mpi
  - 编译时确保启用MPI: cmake -DMPI_FOUND=ON ..
```

# 编译环境快照与回滚

## 功能介绍
该功能可以在编译前自动备份当前环境的关键配置（如环境变量、依赖库版本、CMake配置等），如果编译失败，可以一键回滚到初始状态，避免环境被污染。

## 支持的平台
- Linux/macOS: 使用 `environment_snapshot.sh`
- Windows: 使用 `environment_snapshot.ps1` (PowerShell)

## 使用方法

### Linux/macOS

#### 1. 显示当前环境状态
```bash
./environment_snapshot.sh --action status
```

#### 2. 创建环境快照
```bash
./environment_snapshot.sh --action backup
```
或指定快照文件名：
```bash
./environment_snapshot.sh --action backup --file my_snapshot.txt
```

#### 3. 恢复环境
```bash
./environment_snapshot.sh --action restore
```
或指定快照文件：
```bash
./environment_snapshot.sh --action restore --file my_snapshot.txt
```

### Windows (PowerShell)

#### 1. 显示当前环境状态
```powershell
.\environment_snapshot.ps1 -Action status
```

#### 2. 创建环境快照
```powershell
.\environment_snapshot.ps1 -Action backup
```
或指定快照文件名：
```powershell
.\environment_snapshot.ps1 -Action backup -SnapshotFile my_snapshot.txt
```

#### 3. 恢复环境
```powershell
.\environment_snapshot.ps1 -Action restore
```
或指定快照文件：
```powershell
.\environment_snapshot.ps1 -Action restore -SnapshotFile my_snapshot.txt
```

## 快照内容
快照文件包含以下关键信息：
- 系统信息（操作系统、用户、当前目录）
- 环境变量（CXX_COMPILER, CC_COMPILER, CMAKE_PREFIX_PATH等）
- MPI相关库的安装路径
- CMake版本信息
- 编译器版本信息

## 注意事项
1. 已安装的依赖包无法自动回滚，请手动卸载新增的包
2. 建议在每次编译前创建环境快照
3. 快照文件仅记录当前会话的环境变量，永久环境变量需要手动恢复
4. Windows系统需要使用PowerShell运行脚本
