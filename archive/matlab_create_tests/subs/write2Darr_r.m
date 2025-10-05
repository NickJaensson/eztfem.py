% helper function to write a 2D array of reals to a file
function write2Darr_r(skip,arr,name)
    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines(skip+"["+sprintf('%25.16e,',arr(i,:))+"],");
    end
    mywritelines("    ])");
end
