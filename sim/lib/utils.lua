-- utils.lua
local M = {}

-- Detect OS once
local is_windows = package.config:sub(1,1) == "\\"

-- Internal Command Mapping
M.is_windows = is_windows
M.CP   = is_windows and "copy /y" or "cp"
M.RM   = is_windows and "del /q"  or "rm -f"
M.MOVE = is_windows and "move"    or "mv"

-- 1. Portable Clear
function M.clear_cmd()
    os.execute(is_windows and "cls" or "clear")
end

-- 2. Portable Directory Creation
function M.ensure_dir(dir)
    local cmd = is_windows and 
        string.format('if not exist "%s" mkdir "%s"', dir, dir) or 
        string.format('mkdir -p "%s"', dir)
    os.execute(cmd)
end

-- 3. Portable Bulk Copy
-- Patterns routinely match zero files (e.g. "at*.dat" on tasks that never
-- produce it), so the shell's per-call echo ("N file(s) copied." / "The
-- system cannot find the file specified.") is expected noise, not a
-- diagnostic -- callers don't check the return value. Silence it.
function M.copy_pattern(pattern, destination)
    -- Normalize destination slashes for Windows shells
    local dest = destination
    if is_windows then dest = dest:gsub("/", "\\") end

    local quiet = is_windows and "> nul 2>&1" or "> /dev/null 2>&1"
    local cmd = string.format('%s %s %s %s', M.CP, pattern, dest, quiet)
    return os.execute(cmd)
end

return M