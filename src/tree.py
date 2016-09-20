###############################################################
####### PROCESSING OF TREES ###################################
###############################################################

def addNode(tree):
    id_node = 0
    while tree.has_key(id_node):
        id_node = id_node + 1
    tree[id_node] = ["N"+str(id_node),-1,[],0,"",""]
    return id_node

def getAncestor(tree):
    if tree.has_key("ancestor"):
        return tree["ancestor"]
    else:
        return -1

def setAncestor(tree,node):
    tree["ancestor"] = node
    
def getLength(tree,node):
    return tree[node][3]

def getName(tree,node):
    return tree[node][0]

def setName(tree,node,name):
    tree[node][0] = name

def getSpecies(tree,node):
    return tree[node][5]

def writeSpecies(tree,node,annot):
    tree[node][5] = annot

def getNodes(tree):
    clefs = tree.keys()
    c = 0
    while c < len(clefs):
        if (clefs[c] == "sequence" or
            clefs[c] == "ancestor" or
            len(tree[clefs[c]]) == 0):
            del clefs[c]
        else:
            c = c + 1
    return clefs

def getParent(tree,node):
    return tree[node][1]

def addChild(tree,pere,node):
    tree[pere][2].append(node)

def removeNodeAndChildren(tree,node):
    children = list(getChildren(tree,node))
    for child in children:
        removeNodeAndChildren(tree,child)
    removeNode(tree,node)

def removeNode(tree,node):
#    print "effacement du noeud",node
    tree[node] = []

def removeChildAndChildren(tree,pere,node):
    numero = 0
    while node != getChild(tree,pere,numero):
        numero = numero + 1
    del tree[pere][2][numero]
    removeNodeAndChildren(tree,node)
    
def removeChild(tree,pere,node):
    numero = 0
    while node != getChild(tree,pere,numero):
        numero = numero + 1
    del tree[pere][2][numero]  
    removeNode(tree,node)

def getChild(tree,node,k):
	return tree[node][2][k]

def getNumberOfChildren(tree,node):
    return len(tree[node][2])

def getChildren(tree,node):
    return tree[node][2]

def getBrother(tree,node):
    anc = getParent(tree,node)
    if (getChild(tree,anc,0) == node):
        return getChild(tree,anc,1)
    else:
        return getChild(tree,anc,0)
			
def isLeaf(tree,node):
    return (len(getChildren(tree,node)) == 0)

def isRoot(tree,node):
    return (tree[node][1] == -1)

def isDup(tree,node):
    return (tree[node][4] == "D")

def lastCommonAncestor(tree,a,b):
    ancestor = -1
    ancestorsa = [a]
    while not isRoot(tree,a):
        a = getParent(tree,a)
        ancestorsa.append(a)
    ancestorsb = [b]
    while not isRoot(tree,b):
        b = getParent(tree,b)
        ancestorsb.append(b)
#    print ancestorsa,ancestorsb
    while len(ancestorsa) > 0 and len(ancestorsb) > 0 and ancestorsa[-1] == ancestorsb[-1]:
        ancestor = ancestorsa[-1]
        del ancestorsa[-1]
        del ancestorsb[-1]
#    print "ancestor",ancestor
    return ancestor

def distanceFrom(tree,a,b):
    ancestor = lastCommonAncestor(tree,a,b)
    distance = 0
    while a != ancestor:
        #print tree[a]
        distance = distance + tree[a][3]
        a = getParent(tree,a)
    while b != ancestor:
        #print tree[b]
        distance = distance + tree[b][3]
        b = getParent(tree,b)        
    return distance

def getLeaves(tree,a):
#    print "getleaves",a
    if isLeaf(tree,a):
	return [a]
    else:
        #print "non feuille",child1(a),child2(a)
        result = []
        children = list(getChildren(tree,a))
        for child in children:
            result = result + getLeaves(tree,child)
        return result
    
    
def writeTree(tree,a,NHX):
#    print a,tree[a]
    if isLeaf(tree,a):
        if isRoot(tree,a):
            chaine = "("
        else:
            chaine = ""
	chaine = chaine + tree[a][0] + ":" + str(tree[a][3])
        if tree[a][5] != "":
            chaine = chaine + "[&&NHX:S="+tree[a][5]+"]"
        if isRoot(tree,a):
            chaine = chaine + ")"        
    else:
        chaine = "("
        children = list(getChildren(tree,a))
        for child in children:
            chaine = chaine + writeTree(tree,child,NHX)+","
        chaine = chaine[:-1]+")"
        if not isRoot(tree,a):
            chaine = chaine + ":" + str(tree[a][3])  
        if NHX and tree[a][4] != "" or tree[a][5] != "":
            chaine = chaine + "[&&NHX:"
            if tree[a][5] != "":
                chaine = chaine + "S="+tree[a][5]
            if tree[a][4] == "D" or  tree[a][4] == "WGD":
                chaine = chaine+":D=Y"
            chaine = chaine + "]"
    if isRoot(tree,a):
        chaine = chaine + ";"
    return chaine

def getRoot(tree):
    keys = getNodes(tree)
    start = keys[0]
    while (not isRoot(tree,start)):
        start = getParent(tree,start)
    return start

def getNodesBetween(tree,a,b):
    chemin = []
    ancestor = -1
    ancestorsa = []
    while not isRoot(tree,a):
        a = getParent(tree,a)
        ancestorsa.append(a)
    ancestorsb = []
    while not isRoot(tree,b):
        b = getParent(tree,b)
        ancestorsb.append(b)
    while len(ancestorsa) > 0 and len(ancestorsb) > 0 and ancestorsa[-1] == ancestorsb[-1]:
        ancestor = ancestorsa[-1]
        del ancestorsa[-1]
        del ancestorsb[-1]
#    print "ancestor",ancestor
    return ancestorsa+[ancestor]+ancestorsb


def isAncestor(tree,a,b):
    result = False
    current = b
    while ((not result) and
           (not isRoot(tree,current))):
        if current == a:
            result = True
        else:
            current = getParent(tree,current)
    return result
    
def treeCopy(tree):
	result = {}
	for k in tree.keys():
		if k == "ancestor" or k == "sequence":
			result[k] = tree[k]
		else:
			result[k] = [tree[k][0],tree[k][1],list(tree[k][2]),tree[k][3],tree[k][4],tree[k][5]]
	return result


def changeRoot(tree,newRoot): # the new root is between newRoot and it parent
	if not isRoot(tree,newRoot):
		print "not available for the moment"
	

#####################################################
#####################################################
#  Traversal of one tree 
# 
#####################################################
#####################################################


def readTree(treeseq):

    ###############################################
    ######### TREE READING ########################
    ###############################################
    tree = {"sequence":treeseq}
    id_node = 0
    pile = []
    t = 0
    while t < len(treeseq):
#        print t,treeseq[t],len(pile)
        if treeseq[t] == "(":
            id_node = id_node + 1
            tree[id_node]=["N"+str(id_node),-1,[],-1,"",""]
                        # [nom,pere,[enfants],longueur,annotD,dotannot]
#            print "ouverture",tree[id_node]
            if len(pile) > 0:
                tree[id_node][1] = pile[-1]
            pile.append(id_node)
            t = t + 1
        elif treeseq[t] == ")":
            t = t + 1
            
            if treeseq[t] == "@":
                t = t + 1
                tree["ancestor"] = pile[-1]
    
            while (treeseq[t] != ":" and
                   treeseq[t] != ";" and
                   treeseq[t] != "[" and
                   treeseq[t] != ")" and
                   treeseq[t] != ","):
                t = t + 1
                
                
            if treeseq[t] == ":":
                debut = t + 1
                while treeseq[t] != "," and treeseq[t]!=")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                longueur = float(treeseq[debut:t])
                tree[pile[-1]][3] = longueur
                while treeseq[t] != "," and treeseq[t] != ")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                    
            if treeseq[t] == "[":
                debut = t + 1
                t = debut + treeseq[debut:].find("]")
                chaine = treeseq[debut:t]
                mots = chaine.split(":")
                for m in mots:
                    if m == "D=Y":
                        tree[pile[-1]][4] = "D"
                    if m[:2] == "S=":
                        tree[pile[-1]][5] = m[2:]
                t = t + 1
               
            del pile[-1]
            
            if treeseq[t] == ";":
                t = len(treeseq)
                
        elif treeseq[t] == ";":
            t = len(treeseq)
            
        elif treeseq[t]==",":
            t = t + 1
            
        else:  # nom d'une feuille
            id_node = id_node + 1
            tree[id_node] = ["",-1,[],-1,"",""]
            if len(pile)>0:
                tree[id_node][1]=pile[-1]
            pile.append(id_node)
            debut = t
            while (treeseq[t]!="," and
                   treeseq[t]!=")" and
                   treeseq[t]!=":" and
                   treeseq[t] != "["):
                t=t+1
            nom = treeseq[debut:t]
            tree[pile[-1]][0] = nom
            
            if treeseq[t]==":":
                debut = t + 1
                while treeseq[t]!="," and treeseq[t]!=")" and treeseq[t] != "[" and treeseq[t] != ";":
                    t = t + 1
                longueur = float(treeseq[debut:t])
                tree[id_node][3] = longueur
                
            #print "fin nom"
            if treeseq[t] == "[":
                debut = t + 1
                t = debut + treeseq[debut:].find("]")
                chaine = treeseq[debut:t]
                mots = chaine.split(":")
                for m in mots:
                    if m[:2] == "S=":
                        tree[pile[-1]][5] = m[2:]
                t = t + 1
                
            del pile[-1]

    # remplissage des enfants
    nodes = list(getNodes(tree))
    for node in nodes:
        if not isRoot(tree,node):
	    pere = getParent(tree,node)
            addChild(tree,pere,node)
        
    return tree
