# g-sqz
# supplemental testing helper functions

def count_nucleotides(file_name):
    file = open(file_name, 'r')
    count = 0
    read = True
    while read:
        # sequence identifier
        line_1 = file.readline().rstrip('\n')
        if not line_1:
            read = False
        else:
            start = line_1.index('=')+1
            count += int(line_1[start:len(line_1)])
            line_2 = file.readline()
            line_3 = file.readline()
            line_4 = file.readline()
        file.flush()
    file.close()
    return count

#test data
print(count_nucleotides('test2.fastq'))

