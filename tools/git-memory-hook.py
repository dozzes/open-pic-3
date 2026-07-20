#!/usr/bin/env python3
"""
PostToolUse hook: fired after every Bash tool call.
If the command was a git commit, prints a reminder for Claude to update
the project memory at C:\Users\lamar\.claude\projects\c--Work-open-pic-2\memory\
"""
import json
import sys

try:
    data = json.load(sys.stdin)
    cmd = data.get("tool_input", {}).get("command", "")
    if "git commit" in cmd:
        print(
            "[memory-hook] git commit executed — "
            "check whether MEMORY.md or project memory files need updating."
        )
except Exception:
    pass
