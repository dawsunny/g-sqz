# g-sqz
# The Huffman tree builder class

import _pickle as pickle
from HuffmanNode import *
from struct import *


# reads the file and builds a dictionary of data and its frequency
def build_map(file_name):
    if file_name.endswith('.fastq'):
        huffman_map = build_map_fastq(file_name)
    return huffman_map

# builds a map for FASTQ files
def build_map_fastq(file_name):
    huffman_map = {}
    file = open(file_name, 'r')
    read = True
    while read:
        # sequence identifier
        line_1 = file.readline()
        # checks for an empty file or for the end of a file
        if not line_1:
            read = False
        elif line_1[0] != '@':
            raise FileFormatIncorrectException('The sequence identifier line does not match the FASTQ format')
        else:
            # raw sequence
            line_2 = file.readline().rstrip('\n')
            # line 3 is relatively unimportant given it is a comment (for now)
            file.readline()
            # quality scores
            line_4 = file.readline().rstrip('\n')
            # more of an length error check
            if len(line_2) == len(line_4):
                # count the indices in the line (NOTE: does NOT count '\n')
                for i in range(0, len(line_2)):
                    # store the key as sequence-score
                    key = "" + line_2[i] + line_4[i]
                    if key in huffman_map:
                        huffman_map[key] += 1
                    else:
                        huffman_map[key] = 1
            else:
                raise FileFormatIncorrectException('The length of the raw sequence does not match the length of the quality score')
        file.flush()
    file.close()
    # sorts and returns acc to ascending order of values
    return sorted(huffman_map.items(), key=lambda key: key[1], reverse=True)

# converts the sorted list to Huffman nodes
def generate_huffman_nodes(sorted_list):
    huffman_node_list = []
    for i in sorted_list:
        node = HuffmanNode(None, None, i[0])
        huffman_node_list.append((i[1], node))
    return huffman_node_list

# builds the Huffman tree
def build_huffman_tree(huffman_node_list):
    while (len(huffman_node_list) > 1):
        right_node_tuple = huffman_node_list.pop()
        left_node_tuple = huffman_node_list.pop()
        parent_node = HuffmanNode(left_node_tuple[1], right_node_tuple[1], None)
        insert_node(huffman_node_list, (left_node_tuple[0]+right_node_tuple[0], parent_node))
    return huffman_node_list[0][1]

# inserts the node at the appropriate position
# more effective than sort and checking through all elements (in most cases)
def insert_node(huffman_node_list, node_tuple):
    data = node_tuple[0]
    pointer = len(huffman_node_list)-1
    while (pointer >= 0):
        if data <= huffman_node_list[pointer][0]:
            huffman_node_list.insert(pointer+1, node_tuple)
            return True
        else:
            pointer-=1
    huffman_node_list.insert(0, node_tuple)
    return True     

# TODO everything below this -----

# builds a Huffman Tree dict for easy access of locations
def buildHuffList(filename):
    tree, huffman_map = buildHuffNode(filename)
    # sorts in descending order of freqeuncy - it is IMPORTANT
    huffman_map.sort(key=lambda k: k[1])
    # create lists of nucleotides
    temp = ['G', 'C', 'A', 'T', 'N']
    # builds master list
    hList = [[], [], [], [], []]
    # adds to master list
    # this will append() larger freqs first - optimal access time values!!!
    for i in range(0, len(huffman_map)):
        val = huffman_map[i][0]
        code = findHuffCode(val, tree)
        index = temp.index(val[0])
        hList[index].append([val[1], code])
    return hList, tree


# finds the HuffMan (en)Code for the given value
def findHuffCode(val, tree):
    code = ""
    while (tree.isLeaf() is False):
        data = tree.getLeft().getData()
        if (data.count(val) == 0):  # will give 0 or 1... better than using for loop! =)
            tree = tree.getRight()
            code += "1"
        else:
            tree = tree.getLeft()
            code += "0"
    #print val + " " + code
    return code


def writeHuffFile(filename, newfile):
    hList, tree = buildHuffList(filename)
    # open both files
    new = open(newfile, 'wb')
    old = open(filename, 'r')
    # dump object into new file
    pickle.dump(tree, new)
    nLine = pack('c', '\n')
    new.write(nLine)
    # write the remaining bits
    while True:
        o1 = old.readline()
        if not o1:
            break
        else:
            o2 = old.readline()
            # line 3 is relatively unimportant
            old.readline()
            o4 = old.readline()
            # omits "length = 0" lines ---> pointless!
            if (len(o2) != 1 and len(o4) != 1):
                # read the designated string value
                stt = o1.find(" ") + 1
                stp = o1.find(" ", stt)
                # write the designated string value
                desg = pack('s', o1[stt:stp])
                new.write(desg)
                # write the bits for Huffman
                binCode = ""  # the bitstring
                for i in range(0, len(o2) - 1):
                    nuc = "" + o2[i]
                    val = "" + o4[i]
                    binCode += getBinValue(nuc, val, hList)
                # get and write remainder to file
                rmn = len(binCode) % 8
                rem = chr(rmn)
                new.write(rem)
                # this is epic - makes life a lot easier if DONE HERE
                for i in range(rmn, 8):
                    binCode += "0"
                # intialize x 
                x = 0
                # extracts fragments of 8bits
                for i in range(0, len(binCode)):
                    # extracts indices from fragments
                    c = int(binCode[i])
                    # convert to asci
                    x = x << 1 or c
                    # write bit
                    if (i % 8 == 7):
                        bit = chr(x)
                        new.write(bit)
                        x = 0
                # write new line
                new.write(nLine)
    # close files
    new.close()
    old.close()
    return True


# gets the binary value from the
def getBinValue(n, v, h):
    temp = ['G', 'C', 'A', 'T', 'N']
    subList = h[temp.index(n)]
    for i in range(0, len(subList)):
        if (subList[i][0] == v):
            return subList[i][1]

            # write format:
            # 1. pickled object
            # 2. \n
            # 3. RECURSIVE:
            # a. serial number
            # b. binary remainder
            # c. binary code
            # d. \n

class FileFormatIncorrectException(Exception):
    def __init__(self, error):
        self.error = error
        Exception.__init__(self, 'File Format Incorrect Exception: %s' % error)
