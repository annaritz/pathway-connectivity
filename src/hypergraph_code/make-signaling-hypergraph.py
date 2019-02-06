#/usr/bin/python

import sys
import re
import glob

DELIM=';'


################################
def readFiles(prefix):
    elementfile = '%s-elements.txt' % (prefix)
    reactionfile = '%s-reactions.txt'  % (prefix)
    complexfile = '%s-complexes.txt'  % (prefix)
    controlfile = '%s-controls.txt'  % (prefix)
    subpathwayfile = '%s-subpathways.txt'  % (prefix)
    entitysetfile = '%s-entitysets.txt' % (prefix)
    print 'element file is %s\nreaction file is %s\ncontrol file is %s\n' % (elementfile,reactionfile,controlfile)

    # read elements
    # columns:  #id name altNames altIDs elementtype features cellularLocation evidence
    # indices    0   1     2         3         4            5            6        7
    elementids = readItemSet(elementfile,1) #1-indexed

    ## NEW add entity sets as elements
    elementids.update(readItemSet(entitysetfile,1))

    # read complexes
    # columns: #id name features celularLocation elements stoichiometry evidence
    # indices   0   1      2            3           4          5           6 
    # complexes = {reactomeid: list of reactomeid components}
    complexes = {reactomeid:elements.split(DELIM) for reactomeid,elements in readColumns(complexfile,1,5)} # 1-indexed

    # read pathways
    # columns: #id name xref
    # indices:  0   1    2
    pathwayids = readItemSet(subpathwayfile,1) # 1-indexed
   
    # read controls
    # columns: #id controllingElements controlledElements controlType catalysisdir evidence xref
    # indices:  0          1                    2              3          4           5      6
    lines = [list(l) for l in readColumns(controlfile,1,2,3,4)] # 1-indexed
    controlids = set([n for n,c1,c2,c3 in lines])
    controlling = {n:c1.split(DELIM) for n,c1,c2,c3 in lines}
    controlled = {n:c2.split(DELIM) for n,c1,c2,c3 in lines}
    
    controltypes = {}
    for n,c1,c2,c3 in lines:
        if c3 == "ACTIVATION":
            controltypes[n] = 1
        elif c3 == "ACTIVATION_ALLOSTERIC":
            controltypes[n] = 1
        elif c3 == "INHIBITION":
            controltypes[n] = 0
        elif c3 == "INHIBITION_COMPETITIVE":
            controltypes[n] = 0
        elif c3 == "INHIBITION_NONCOMPETITIVE":
            controltypes[n] = 0
        elif c3 == "INHIBITION_ALLOSTERIC":
            controltypes[n] = 0
        else:
            sys.exit('ERROR: control type %s not specified\n' % (c3))
    controlledelements = set()
    for n in controlled:
        controlledelements.update(controlled[n])

    # read reactions
    # columns: #id E F R interactionType E&Fstoichiometry spontaneous conversiondirection evidence xref
    # indices:  0  1 2 3       4              5               6                7             8      9
    lines = [list(l) for l in readColumns(reactionfile,1,2,3,4)] #1-indexed
    reactionids = set([n for n,e,f,r in lines])
    tails = {n:e.split(DELIM) for n,e,f,r in lines}
    heads = {n:f.split(DELIM) for n,e,f,r in lines}
    regs = {n:r.split(DELIM) for n,e,f,r in lines}

    ## CHECK 0: Remove 'TemplateReactions' from tail.
    for rid in reactionids:
        if 'TemplateReaction' in rid:
            tails[rid] = ['None']
             
    ### CHECK 1: Make sure that all controlled elements have a regulator.
    for n in controlledelements:
        ## get elements that CONTROL this element
        controlids = [cid for cid in controlled if n in controlled[cid]]
        ## to flatten: [item for sublist in l for item in sublist]
        controllingids = [controlling[cid] for cid in controlids]
        controllingids = [item for sublist in controllingids for item in sublist]

        if n in regs and set(controlids) == set(regs[n]): ## we're good! continue
            continue

        if 'TemplateReaction' in n: # this is a template reaction that's regulated
            ## n is also a reactionid
            if n not in reactionids:
                sys.exit('ERROR: template reaction %s should be a reaction id' % (n))
            if regs[n] != ['None']:
                sys.exit('ERROR: template reaction %s has a regulator already!' % (n))
            regs[n] = controlids
        elif 'Pathway' in n: # this is a pathway! make new reaction that is directly regulated.
            for c in controlids:
                rid = 'new_'+c
                reactionids.add(rid)
                tails[rid] = ['None']
                heads[rid] = [n]
                regs[rid] = controlling[c]
                controltypes[rid]=controltypes[c]
        elif 'Catalysis' in n:
            for c in controlids:
                if 'Modulation' not in c:
                    sys.exit('ERROR: controlling element %s is not listed as a regulator; but controlling ID %s is not a modulator.' % (n,c))
                # this controls a control!! Make a new reaction that directly regulates.
                # the controller of c (controlling[c]) controls the controller of n (controlling[n])
                rid = 'new_'+c
                reactionids.add(rid)
                tails[rid] = ['None']
                ## set the controlling element of the control element (n)  as the head.
                heads[rid] = controlling[n]
                regs[rid] = [c]
                controltypes[rid]=controltypes[c]
                controlled[c].append(rid)
        else:
            print 'n=%s' % (n)
            #print 'controlids=',controlids
            #print 'controllingids=',controllingids
            #print 'tails[n]=',tails[n]
            #print 'heads[n]=',heads[n]
            #print 'regs[n]=',regs[n]
            print 'UNACCOUNTED CHECK!!'
            #sys.exit()
            ## currently 124 unaccounted checks! Continue.
            
    ### CHECK: it's possible that a reaction has no T,R+,or R-.  
    ## For example, if it's a template reaction and there's no control element for it.  
    ## Remove these cases.
    toremove = set()
    for n in reactionids:
        if (tails[n] == ['None'] and regs[n] == ['None']):
            print 'WARNING: %s has an empty tail and regulator set.' % (n)
            toremove.add(n)
    for n in toremove:
        reactionids.remove(n)
        del tails[n]
        del heads[n]
        del regs[n]

    ### CHECK 2: Make sure all controls are present in reactions. (used for NCI-PID)
    ### e.g. when A inhibits B, this will be in controls but not necessarily in reactions.
    return elementids,pathwayids,complexes,controlling,controlled,controltypes,reactionids,tails,heads,regs
        
################################
def writeFiles(elementids,pathwayids,complexes,controlling,controlled,controltypes,reactionids,tails,heads,regs,prefix):
    eout = open('%s-hyperedges.txt' % (prefix),'w')
    eout.write('#Tail\tHead\tPosReg\tNegReg\tID\n')

    numskipped = 0
    numhyperedges = 0
    for n in reactionids: # go through all reaction REACTOMEIDs
        if heads[n] == ['None']:
            print 'WARNING: reaction %s has an empty head. Skipping.' % (n)
            numskipped+=1
            continue
    
        ## TODO: check for self loops?

        ## divide regulators into positive regs and negative regs
        posregs = set()
        negregs = set()
        if regs[n] != ['None']:

            if heads[n][0] in pathwayids:
                print 'Writing directly-regulating edge (%s in pathwayids)' % (heads[n][0])
                if controltypes[n] == 1:
                    posregs.update([r for r in regs[n]]) 
                else:
                    negregs.update([r for r in regs[n]])                     
            
            else:
                for r in regs[n]:
                    if r in controlled and n in controlled[r]:
                        ## NEW: there are some control/catalysis that are controlled by None? E.g. Control289 in Reactome.
                        if r in controltypes and controltypes[r] == 1: # ACTIVATION
                            posregs.update([a for a in controlling[r] if a != 'None']) ## Add the elements that CONTROL the regulation!
                        elif r in controltypes and controltypes[r] == 0: # INHIBITION
                            negregs.update([a for a in controlling[r] if a != 'None']) ## Add the elements that CONTROL the regulation!
                        else:
                            ## TODO; originally, if no controltypes[r] then put ACTIVATION
                            print 'ERROR: r not in controltypes!'
                            print 'r=',r
                            print 'n=',n
                            sys.exit()
                    else:
                        ## TODO: originally if no controlled[r] then put 'UnkownMechanism'.
                        print 'n=',n
                        print 'tails[n]=',tails[n]
                        print 'heads[n]=',heads[n]
                        print 'regs[n]=',regs[n]
                        print 'ERROR: r not in controlled or n not in controlled[r]'
                        print 'r=',r
                        print 'r in controlled:',r in controlled
                        if r in controlled:
                            print controlled[r]
                            print 'n in controlled[r]:',n in controlled[r]
                        sys.exit()

        # len of these sets may be 0 of regs[n] == ['None']
        # or there are either all pos or all neg regulators
        # (the other list will be empty).
        if len(posregs)==0:
            posregs.add('None')
        if len(negregs)==0:
            negregs.add('None')

        eout.write('%s\t%s\t%s\t%s\t%s\n' % (DELIM.join(tails[n]),DELIM.join(heads[n]),DELIM.join(posregs),DELIM.join(negregs),n))
        numhyperedges+=1
    
    ## TODO: if a complex has 0 incoming edges, should we make the proteins its incoming edges?

    eout.close()
     
    nout = open('%s-hypernodes.txt' % (prefix),'w')
    nout.write('#hypernode\tnodes\n')
    numhypernodes = 0

    torecheck = set() # these are complexes that contain other complex IDs.

    for n in complexes: # write complex IDs as multiple elements.
        ## check that all elements in the complex are in elementids
        ## if some member is ANOTHER complex, "flatten" it to the base elements.
        elementsincomplex = set()
        for c in complexes[n]:
            elementsincomplex.update(getComplexMembers(c,complexes,elementids,verbose=False))

        #print 'writing %s' % (n)
        nout.write('%s\t%s\n' % (n,DELIM.join(elementsincomplex)))
        numhypernodes+=1

    ## we write element IDs AFTER hypernodes.
    ## NOTE THIS also write EntitySets as singleton hypernodes
    for n in elementids: # write element IDs as singleton hypernodes
        nout.write('%s\t%s\n' % (n,n))
        numhypernodes+=1

    for n in pathwayids: # write pathwayIDs as singleton hypernodes
        nout.write('%s\t%s\n' % (n,n))
        numhypernodes+=1
 
    nout.close()

    print '\n'
    print '%d reactions written' % (numhyperedges)
    print '%d reactions skipped' % (numskipped)
    print '%d hypernodes written (%d nodes, %d complexes, %d subpathways)' % (numhypernodes,len(elementids),len(complexes),len(pathwayids))

    print '\n'
    print 'Wrote to %s-hyperedges.txt' % (prefix)
    print 'Wrote to %s-hypernodes.txt' % (prefix)
    return

################################
def getComplexMembers(el,complexes,elementids,verbose=False):
    if verbose:
        print 'Checking Element',el
    currentset = set()
    if el == 'None':
        return currentset
    if el in elementids:
        currentset.add(el)
    elif el in complexes: ## It's a complex!
        for c in complexes[el]:
            currentset.update(getComplexMembers(c,complexes,elementids,verbose))
    else:
        print 'ERROR: %s is not a complex or an elementID' % (el)
        sys.exit()
    return currentset
 
################################  
## from utilsPoirel.py

def readItemList(f, col=1, sep='\t'):
    '''
    Read the given column of the tab-delimited file f
    and return it as a list. Col is the 1-based column index.
    '''
    itemlist = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split(sep)
        if len(items)<col:
            continue
        itemlist.append(items[col-1])
    return itemlist
    
def readItemSet(f, col=1, sep='\t'):
    '''
    Read the given column of the tab-delimited file f
    and return it as a set. Col is the 1-based column index.

    A wrapper to readItemList, returning a set instead of a list.
    '''
    return set(readItemList(f, col=col, sep=sep))

def readColumns(f, *cols):
    '''
    Read multiple columns and return the items from those columns
    in each line as a tuple.

    foo.txt:
        a b c
        d e f
        g h i

    Calling "readColumns('foo.txt', 1, 3)" will return:
        [(a, c), (d, f), (g, i)]

    '''
    if len(cols)==0:
        return []
    rows = []
    for line in open(f, 'r').readlines():
        if line=='':
            continue
        if line[0]=='#':
            continue
        items = line.rstrip().split('\t')
        if len(items)<max(cols):
            continue
        rows.append(tuple([items[c-1] for c in cols]))
    return rows
                  
################################
def main(args):

    if len(args) != 3:
         print 'USAGE: make-signaling-hypergraph.py <DATADIRECTORY> <OUTDIR>\n\t<DATADIRECTORY>: directory containing the files parsed by BioPAX Paxtools parser called ReactomeParser.jar\n\t<OUTDIR>: output directory for signaling hyergraph files.\n\n'
         sys.exit()

    DATADIR = args[1]
    DIR = args[2]

    files = glob.glob(DATADIR+'*-reactions.txt')
    print '%d files' % (len(files))
    c=1
    for f in files:
        dataprefix = f[:-14]
        print 'FILE %d:' % (c),dataprefix
        c+=1
        elementids,pathwayids,complexes,controlling,controlled,controltypes,reactionids,tails,heads,regs = readFiles(dataprefix)

        outprefix = DIR+dataprefix.split('/')[-1]
        writeFiles(elementids,pathwayids,complexes,controlling,controlled,controltypes,reactionids,tails,heads,regs,outprefix)

if __name__=='__main__':
    main(sys.argv)
