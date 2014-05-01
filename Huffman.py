import _pickle as pickle
from HuffNode import *
from struct import *

# builds a list of data and its frequency
def buildMap(filename):
    hMap = {}
    f = open(filename, 'r')
    while True:
        l1 = f.readline()
        if not l1:
            break
        else:
            l2 = f.readline()
            # line 3 is relatively unimportant
            f.readline()
            l4 = f.readline()
        if len(l2) == len(l4):
            # count the indices in the line (NOTE: does NOT count '\n')
            for i in range(0, len(l2) - 1):
                key = "" + l2[i] + l4[i]
                if hMap.has_key(key):
                    val = hMap[key]
                    val += 1
                else:
                    val = 1
                hMap[key] = val
        f.flush()
    f.close()
    # sorts and returns acc to ascending order of values
    return sorted(hMap.items(), key=lambda k: k[1])


# builds the Huffman Tree
def buildHuffNode(filename):
    # build the file
    hMap = buildMap(filename)
    # create node map
    nMap = []
    # create a list of HuffMan Tree from file
    for i in range(0, len(hMap)):
        val = hMap[i]
        node = HuffNode(None, None, val[0], val[1])
        # add to nMap
        nMap.append([val[1], node])
    while True:
        if len(nMap) == 1:
            break
        else:
            l = nMap.pop(0)
            r = nMap.pop(0)
            par = HuffNode(l[1], r[1], None, None)
            nMap.append([par.getFreq(), par])
            # sorts the map with new values
            nMap.sort()
    t = nMap.pop(0)
    tree = t[1]
    return tree, hMap


# builds a Huffman Tree dict for easy access of locations
def buildHuffList(filename):
    tree, hMap = buildHuffNode(filename)
    # sorts in descending order of freqeuncy - it is IMPORTANT
    hMap.sort(key=lambda k: k[1])
    # create lists of nucleotides
    temp = ['G', 'C', 'A', 'T', 'N']
    # builds master list
    hList = [[], [], [], [], []]
    # adds to master list
    # this will append() larger freqs first - optimal access time values!!!
    for i in range(0, len(hMap)):
        val = hMap[i][0]
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
