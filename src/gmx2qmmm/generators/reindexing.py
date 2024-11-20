def reindexing_memory(ref, outeratomlist): # Added by Nicola
    c = 0
    outeratomlist = set([int(i) for i in outeratomlist])
    mem = []
    for j in range(len(ref) - 1):
        if int(ref[j]) in outeratomlist:
            mem.append([ref[j + 1], ref[j]])
            if mem[-1][0] != mem[-1][1]:
                c += 1
                mem[-1] = [mem[-1][1], mem[-1][0] - c]
            else:
                continue
        else:
            if not mem == [] and mem[-1][0] != mem[-1][1]:
                mem.append([ref[j + 1], ref[j]])
                if mem[-1][0] != mem[-1][1]:
                    mem[-1] = [mem[-1][1], mem[-1][0] - (c + 1)]
        
            else:
                mem.append([ref[j], ref[j]])
                continue


    if int(ref[-1]) not in outeratomlist:
        if int(ref[-2]) not in outeratomlist:
            mem.append([int(ref[-1]), mem[-1][1] + 1])
        else:
            mem.append([int(ref[-1]), mem[-1][1]])
    else:
        mem.append([int(ref[-1]), mem[-1][1]])

    memory = [] #[[old index, new index]]
    for i in mem:
        if int(i[0]) not in outeratomlist:
            memory.append(i)
        else:
            continue
    c = 1
    for i in range(len(memory)):
        memory[i][1] = c
        c += 1 
              
    memory_dict = {item[0]: item[1] for item in memory}

    return memory_dict


def reindexing_memory2(ref, outeratomlist):
    if not ref:
        return {}

    c = 0
    outeratomlist = set(int(i) for i in outeratomlist)
    mem = []

    for j in range(len(ref) - 1):
        current = int(ref[j])
        next_value = int(ref[j + 1])

        if current in outeratomlist:
            mem.append([next_value, current])
            if mem[-1][0] != mem[-1][1]:
                c += 1
                mem[-1][1] -= c
        else:
            if mem and mem[-1][0] != mem[-1][1]:
                mem.append([next_value, current])
                if mem[-1][0] != mem[-1][1]:
                    mem[-1][1] -= (c + 1)
            else:
                mem.append([current, current])

    # Handle last element
    last = int(ref[-1])
    if last not in outeratomlist:
        if int(ref[-2]) not in outeratomlist:
            mem.append([last, mem[-1][1] + 1])
        else:
            mem.append([last, mem[-1][1]])
    else:
        mem.append([last, mem[-1][1]])

    # Create memory dictionary
    memory = []
    for i in mem:
        if i[0] not in outeratomlist:
            memory.append(i)

    memory_dict = {item[0]: item[1] for item in memory}
    return memory_dict

