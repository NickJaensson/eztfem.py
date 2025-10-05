% helper function to write a 2D array of integers to a file
function write2Darr_i(skip,arr,name)
    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines(skip+"["+sprintf('%12i,',arr(i,:))+"],");
    end
    mywritelines("    ],dtype=int)");
end