class HuffNode:
    # This data will exist in all
    Name = "HuffNode"
    # __init__ is a class constructor
    # __****__ is usually a special class method.
    def __init__(self, l, r, d, f):
        # These values are created
        # when the class is instantiated.
        if (d is None and f is None):
            self.data = "" + l.data + r.data
            self.freq = 0
            self.left = l
            self.right = r
            self.leaf = False
        elif (l is None and r is None):
            self.data = d
            self.freq = f
            self.left = None
            self.right = None
            self.leaf = True

    def isLeaf(self):
        return self.leaf

    def getLeft(self):
        return self.left

    def getRight(self):
        return self.right

    def getData(self):
        return self.data        
