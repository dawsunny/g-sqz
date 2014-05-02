# g-sqz
# The node of a Huffman tree

class HuffmanNode:
    # class constructor
    # better to have a list of data than a concatenated string in case on DNA represented numerically
    def __init__(self, left, right, data):
        # leaf node
        if (left is None):
            self.left = None
            self.right = None
            self.data = [data]
        # parent node
        else:
            self.left = left
            self.right = right
            self.data = []
            self.data.extend(self.left.data)
            self.data.extend(self.right.data)

    # basic getters
    def is_leaf(self):
        return (self.left == None and self.right == None)

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def get_data(self):
        return self.data        
