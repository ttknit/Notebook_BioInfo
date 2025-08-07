#!/bin/bash
#
# 这个脚本负责为所有任务添加统一的日志和trap机制

# 定义日志文件路径
LOG_FILE="/path/to/your/auto_submit.log"
TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")

# 检查是否提供了要运行的脚本
if [ -z "$1" ]; then
    echo "$TIMESTAMP - ERROR: Wrapper script called without a task script." >> "$LOG_FILE"
    exit 1
fi

# 获取原始任务脚本的名称
TASK_SCRIPT_NAME=$(basename "$1")
TASK_JOB_ID="$SLURM_JOB_ID"

# 定义 trap 函数，用于在脚本退出时记录状态
function log_exit_status {
    local exit_status=$?
    if [ $exit_status -eq 0 ]; then
        echo "$TIMESTAMP - TASK_STATUS: $TASK_SCRIPT_NAME($TASK_JOB_ID) - SUCCESS" >> "$LOG_FILE"
    else
        echo "$TIMESTAMP - TASK_STATUS: $TASK_SCRIPT_NAME($TASK_JOB_ID) - FAILED (exit code $exit_status)" >> "$LOG_FILE"
    fi
}
# 使用 trap 捕获 EXIT 信号
trap log_exit_status EXIT

echo "$TIMESTAMP - STARTING: $TASK_SCRIPT_NAME($TASK_JOB_ID)" >> "$LOG_FILE"

# --- 执行原始任务脚本 ---
# 使用 "$@" 传递所有参数，包括原始脚本路径和它可能需要的任何参数
# 这里使用 source，因为它允许被包装脚本继承当前环境
source "$@"

# 脚本执行到这里时，trap会自动被触发
