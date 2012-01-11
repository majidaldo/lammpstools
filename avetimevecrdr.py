





class vecrdr(object):
    def __init__(self,filepath):
        self.gethdrs(filepath)
        self.f=file(filepath)
        return
    
    def gethdrs(self,fp):
        f=file(fp)
        self.title=f.next()[2:-1]
        f.next()
        self.cols=f.next()[2:-1].split(' ')
        f.close()
        return
    
    
    def gotonextsec(self):
        for aline in self.f:
            sl=aline[:-1].split(' ') #-1 to get rid of endofline char
            if len(sl)==2 and sl[0].isdigit()==True and sl[1].isdigit()==True:
                return sl
    
    def getnextchunk(self):
        try: ts,nrows=self.gotonextsec() #outputs a none at eof
        except: raise StopIteration
        nrows=int(nrows);ts=int(ts)
        chunk=[self.f.next() for aline in xrange(nrows)]
        return ts, chunk
        
    def chunk2vecs(self,chunk):
        """input a section of vector out avetime..out vecs"""
        sls=[str.split(aline[:-1],' ') for aline in chunk]
        return zip(*sls)[1:]
    
    def getnextvecs(self):
        ts,chunk=self.getnextchunk()
        return ts,self.chunk2vecs(chunk)
    
    def getnextnvecs(self,n):
        tss=[];vss=[]
        i=0
        while i<n:
            try:
                ts,vs=self.getnextvecs()
                tss.append(ts);vss.append(vs)
                i+=1
            except StopIteration:
                break
        return tss,zip(*vss) #so that 1st indx is the vec, 2nd one is the sequence
    
    def getnextvecsi(self):
        while True:
            yield self.getnextvecs() 
            #for some reason, have to use gen expr 

import itertools
class vecfilesrdr(vecrdr):
    def __init__(self,filelist):
        if type(filelist)==str:filelist=[filelist]
        filelist.sort()
        super(vecfilesrdr,self).__init__(filelist[0])
        
        #replace single file iter     
        fobjs=itertools.imap(open,filelist)
        self.f=itertools.chain.from_iterable(fobjs)
        return
        
#def returnchunks(filepath):
#    f=file(filepath)
#    for aline in f:
#        if
    
    
    