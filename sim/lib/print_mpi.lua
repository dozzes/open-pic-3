--[[
NAME
  print_mpi

FUNCTIONS
  print_root
  fprint_root
    
NOTES
    Package for printing to standard output (stdout) and files.
    Updated for Lua 5.4 compatibility.
]]

local P = {}

-- In modern Lua, we typically just return the table at the end, 
-- but keeping your global assignment logic:
if _REQUIREDNAME == nil then
    print_mpi = P
else
    _G[_REQUIREDNAME] = P
end

-- Use '...' directly and table.pack to handle arguments in Lua 5.x
-- No implicit newline: callers pass explicit "\n". An automatic one here
-- doubled every line break and produced blank-line-riddled startup output.
function P.print_root(proc_idx, ...)
    if proc_idx == 0 then
        local args = table.pack(...)
        for i = 1, args.n do
            io.write(tostring(args[i]))
        end
    end
end
   
function P.fprint_root(proc_idx, out_file, ...)
    if proc_idx == 0 and out_file then 
        local args = table.pack(...)
        for i = 1, args.n do
           out_file:write(tostring(args[i]))
        end
    end
end

return P