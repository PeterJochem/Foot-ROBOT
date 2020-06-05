function find=find_index_extract(start_index,terminal_condition,searching_list)
    
    while searching_list(start_index)>terminal_condition
        if (start_index)==length(searching_list)
            break;
        end
        start_index=start_index+1;
    end
    find=start_index;
end