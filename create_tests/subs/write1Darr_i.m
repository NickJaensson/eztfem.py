% helper function to write a 1D array of integers to a file
function write1Darr_i(skip,arr,name)
    mywritelines("    "+name+" = np.array([");
    for i=1:length(arr)
        mywritelines(skip+sprintf('%12i,',arr(i)));
    end
    mywritelines("    ],dtype=int)");
end