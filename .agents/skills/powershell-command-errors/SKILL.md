---
description: common command execution errors and gotchas in Windows PowerShell for Antigravity
---

# Command Execution Pitfalls & Skills (Windows PowerShell)

As an AI agent operating in a Windows PowerShell environment, you will often need to execute commands using the `run_command` tool. During previous program executions (such as spawning parallel Julia sweeps in `run_all_sweeps.ps1` and `run_sweeps.bat`), several command-related errors and stalls occurred. 

This document summarizes those errors and provides strict rules to prevent them in the future.

## 1. Backgrounding Tasks and I/O Stalls
**The Error**: Running multiple background processes using `Start-Process -NoNewWindow` causes their standard output and error streams to interleave in the same hidden console. This often leads to I/O blocking, causing the processes (like Julia parametric sweeps) to completely stall.
**The Fix**: 
- If you must use PowerShell to spawn background processes and you don't need to see the window, explicitly redirect streams to files:
  `Start-Process -FilePath "julia" -ArgumentList ... -RedirectStandardOutput "out.log" -RedirectStandardError "err.log" -WindowStyle Hidden`
- Alternatively, use the `run_command` tool's native backgrounding capability by using a small `WaitMsBeforeAsync` instead of relying on PowerShell's `Start-Process`, and capture the output via the `command_status` tool.

## 2. Directory Navigation (`cd` command)
**The Error**: Attempting to use `cd C:\path` using the `run_command` tool. Because each tool call is stateless, the directory change does not persist to the next tool call.
**The Fix**: **NEVER PROPOSE A `cd` COMMAND**. Always use the `Cwd` parameter in the `run_command` tool to set the working directory for that specific command.

## 3. Execution Policy for `.ps1` Scripts
**The Error**: Running a PowerShell script directly (e.g., `.\script.ps1`) fails due to Windows default Execution Policies restricting script execution.
**The Fix**: When you need to execute a `.ps1` script, bypass the execution policy for that specific process:
`powershell.exe -ExecutionPolicy Bypass -File .\script.ps1`

## 4. Default Output Encoding (UTF-16 vs UTF-8)
**The Error**: In Windows PowerShell 5.1, using the `>` operator to redirect output to a file creates a UTF-16LE encoded file by default. This causes parsing errors for downstream scripts (like Python or Julia) that expect UTF-8.
**The Fix**: Use `Out-File` with the `-Encoding utf8` flag instead of the `>` operator, or set the console encoding before running commands that produce native text files:
`command | Out-File -FilePath out.txt -Encoding utf8`

## 5. Linux Utilities vs PowerShell Cmdlets
**The Error**: Using Linux commands (`ls`, `rm`, `cat`, `grep`) within PowerShell command strings. While PowerShell aliases some of these (like `cat` -> `Get-Content`, `ls` -> `Get-ChildItem`), their flags and output formats are entirely different, leading to syntax errors or unexpected array objects instead of plain text.
**The Fix**: 
- For file manipulation/searching, prioritize Antigravity native tools (`list_dir`, `view_file`, `grep_search`, `read_url_content`, `replace_file_content`).
- If you absolutely must use shell commands, use the explicit PowerShell cmdlets (e.g., `Select-String`, `Remove-Item`) and be aware that they return objects, not raw strings. 
- **CRITICAL**: Never run `cat` inside a bash/powershell command to create a new file or append to an existing file. Use `write_to_file`.

## 6. String Interpolation (Quotes)
**The Error**: Using single quotes (`'`) when trying to interpolate environment variables or expressions. In PowerShell, single quotes are literal, and double quotes (`"`) expand variables.
**The Fix**: Use double quotes when passing arguments that need expansion (e.g., `"$ENV:USERPROFILE\Documents"`), but watch out for escaping rules if passing JSON or embedded scripts.
