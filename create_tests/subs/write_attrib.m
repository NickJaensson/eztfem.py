% helper function to write attributes of a struct to a file
function write_attrib(skip,struct,name)
    fns = fieldnames(struct);
    for k=1:numel(fns)
        tmp = struct.(fns{k});
        if isnumeric(tmp)
            if ndims(tmp) > 3
                error("error write_attrib: dims("+fns{k}+") > 3")
            end
            if isscalar(tmp)
                mywritelines(skip+name+"."+fns{k}+" = "+string(tmp));
            elseif all(mod(tmp,1)<1e-15,'all') % array of integers
                if ndims(tmp) == 2 
                    if size(tmp,1)==1 || size(tmp,2)==1
                        write1Darr_i(skip,tmp,name+"."+fns{k});
                    else 
                        write2Darr_i(skip,tmp,name+"."+fns{k});
                    end
                else
                    write3Darr_i(skip,tmp,name+"."+fns{k});
                end
            else 
                if ndims(tmp) == 2 
                    if size(tmp,1)==1 || size(tmp,2)==1
                        write1Darr_r(skip,tmp,name+"."+fns{k}); 
                    else
                        write2Darr_r(skip,tmp,name+"."+fns{k});
                    end
                else
                    write3Darr_r(skip,tmp,name+"."+fns{k});
                end
            end
        end
    end
end