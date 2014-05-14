# g-sqz
# supplemental testing helper functions

from HuffmanTree import *

def theoretical_size(file_name):
    print(file_name)
    huffman_map, line_len = build_map(file_name)
    huffman_node = build_huffman_tree(huffman_map)
    huffman_encode_map = generate_huffman_code_map(huffman_node)
    pickled_map = dumps(huffman_encode_map)
    print('Pickled Map Size: ' + str(len(pickled_map))+'b')
    count = 0
    for key in huffman_map:
        count += huffman_map[key]*len(huffman_encode_map[key])
    count = round(count/8)
    print('Huffman Size: ' + str(count)+'b')
    total = len(pickled_map) + count + 7
    if total > 1048576:
        print('Total Size: ' + str(round(total/1048576))+'MB')
    elif total > 1024:
        print('Total Size: ' + str(round(total/1024))+'KB')
    else:
        print('Total Size: ' + str(total)+'B')
    print()


#test data
theoretical_size('./test/test0.0.fastq')
theoretical_size('./test/test0.1.fastq')
