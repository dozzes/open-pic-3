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
function M.copy_pattern(pattern, destination)
    -- Normalize destination slashes for Windows shells
    local dest = destination
    if is_windows then dest = dest:gsub("/", "\\") end
    
    local cmd = string.format('%s %s %s', M.CP, pattern, dest)
    return os.execute(cmd)
end

return M