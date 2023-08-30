% helper function to write a 1D array of reals to a file
function write1Darr_r(skip,arr,name)
    mywritelines("    "+name+" = np.array([");
    for i=1:length(arr)
        mywritelines(skip+sprintf('%25.16e',arr(i))+",");
    end
    mywritelines("    ])");
end