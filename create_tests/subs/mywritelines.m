% helper function to write some lines of text to a file
function mywritelines(str)
    global fn
    writelines(str,fn,WriteMode="append");
end