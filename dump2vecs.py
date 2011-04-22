
##Author: Majid al-Dosari
#Copyright (c) 2011, Majid al-Dosari
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the <organization> nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



import itertools
import collections
import operator

class dump2vecs(object):
    
    def __init__(self,dumpfile):
        self.dumpfile=open(dumpfile,'r')
        dim=self.getdim();self.dumpfile.seek(0)
        self.box={'x':None,'y':None};
        if dim==3: self.box.update({'z':None})
        self.broken=None# a place for when the natoms in a frame changes
        #self.tsi=[]
        self.ts=None
        self.attribs=self.getatomattribs();self.dumpfile.seek(0)
        #self.aids=self.getatomids1();self.dumpfile.seek(0)
        return
        
    def getdim(self):
        self.scantoboxbounds().next()
        d=0
        for al in self.dumpfile:
            if 'ITEM' in al: break
            else:d+=1
        return d
    
    def getatomattribs(self):
        headings=str.split(self.scantonexttable().next())
        attribs=dict(zip(headings,range(-2,len(headings)-2))) #index to a line
        attribs.pop('ITEM:');attribs.pop('ATOMS')
        #attribs.pop('id')
        return attribs
        
#    def getatomids1(self):
#        fo=self.dumpfile
#        cont=True
#        atomids=[]
#        #get no. of atoms
#        while cont==True:
#            al=fo.readline()
#            if 'ITEM: NUMBER OF ATOMS' in al:
#                #read the no. in the next line
#                natoms=int(fo.readline())        
#                cont=False
#        if natoms==0: return atomids
#        #then read to titles
#        cont=True
#        while cont==True:
#            al=fo.readline()
#            if 'ITEM: ATOMS' in al: cont=False
#        #so the next line should have an atom
#        cont=True;atomcount=0
#        while atomcount<natoms:
#            try:al=fo.readline()
#            except EOFError: raise KeyError, '1st frame incomplete'
#            try:anatomid=al.split()[0]
#            except: raise Exception, '1st frame incomplete'
#            atomcount+=1;atomids.append(int(anatomid))
#        #fo.close()
#        return atomids
    
#    def scantonexttable2(self):
##        for aline in self.dumpfile:
##            if 'ITEM: ATOMS' in aline: return self.dumpfile.tell()
#        while 'ITEM: ATOMS' not in self.dumpfile.readline(): continue
#        return self.dumpfile.tell()
    def scantonexttable(self): return self.scantostr('ITEM: A')#TOMS')
    #faster than scatohead idk why!!    
    
#    def scantonextframe2(self):
##        for aline in self.dumpfile:
##            if 'ITEM: ATOMS' in aline: return self.dumpfile.tell()
#        while 'ITEM: TIMESTEP' not in self.dumpfile.readline(): continue
#        return self.dumpfile.tell()
    def scantonextframe(self): return self.scantostr('ITEM: T')#IMESTEP')
    
    def scantonatoms(self):return self.scantostr('ITEM: N')#ER OF ATOMS
    def scantoboxbounds(self):return self.scantostr('ITEM: B')#OX BOUNDS

    


#if atoms not found and we got to the next frame return none
#self.jump=False; 

#make a C loop that filters
#[ aline if "atmid ' in it]
#then map to extract
#[str.split(al)[3] for al in d2vt.dumpfile if '10 ' in al]
    def scantostr(self,astr):#faster. idk why
        for aline in self.dumpfile:
            if astr not in aline: continue #how is this faster?!?!
            else:#it found the str
#                if astr not in self.jumperi: self. #deque
#                else: self.jumperi
                yield aline #self.jump=False
    def scantohead(self,astr):# faster yet
        nc=len(astr)
        for aline in self.dumpfile:
            if astr in aline[:nc]: yield aline

    #so enable jump before each
        #while astr not in self.dumpfile.readline(): continue 
        #return self.dumpfile.tell() 
    
    def fastscantonextatomline(self,atomid): #this is faster
    #than putting it in a list compr
        """atom id has to be first thing in line"""
        aidstr=str(atomid)+' '        
        for ats in self.scantonextframe():
            yield self.scantohead(aidstr).next()
        
            
#slow            
#            try:
#                aid=(int(al.split()[0])) #todo AND
#                if aid== atomid: yield al
#            except:pass
            
        #return self.scantostr(str(atomid)+' ')
        
    def gotoframe(self,ts):#useless? already visited
        self.dumpfile.seek(self.tsi[ts]); return
    
    def gotopreviousframe(self):pass

    def nexttimestep(self):
#        for ani in self.scantonextframe(): #todo. make this jump
#            self.curts=int(self.dumpfile.next())
#            yield self.curts
            #dont yield int(theresult) , it's ok here b/c it wont spend 
            #much time here
        self.scantonextframe().next()
        #self.tsi.append(int(self.dumpfile.next()))
        self.ts=int(self.dumpfile.next())
        return self.ts#[:-1]
    def nextnatoms(self):
#        for ani in self.scantonatoms():
#            self.natoms=int(self.dumpfile.next())
#            yield self.natoms
        self.scantonatoms().next()
        self.natoms=int(self.dumpfile.next())
        return self.natoms
    
    def nextboxbounds(self):
        def txt2tpl(txt): return tuple(map(float,tuple(txt.split())))
#        for ani in self.scantoboxbounds():
#            self.box['x']=txt2tpl(self.dumpfile.next())
#            self.box['y']=txt2tpl(self.dumpfile.next())
#            if len(self.box)==3:
#                self.box.update({'z':txt2tpl(self.dumpfile.next())})
#        yield self.box
        self.scantoboxbounds().next()
        self.box['x']=txt2tpl(self.dumpfile.next())
        self.box['y']=txt2tpl(self.dumpfile.next())
        if len(self.box)==3:
            self.box.update({'z':txt2tpl(self.dumpfile.next())})
        return self.box
        
    def nextframeinfoi(self):
        while True:
            self.nexttimestep()
            self.nextnatoms()
            self.nextboxbounds()#.next()
            self.scantonexttable().next()
            yield self.ts, self.box.values() ,self.natoms \
            , iter(([self.dumpfile.next() for al in xrange(self.natoms)])) #atom lines
        #so this errors out if you're at the end and it's incomplete        
        
        #what happens when lines are less AND at EoF?
        #A:you get StopIter meaning the frame is incomplete
        #so the frame should be completed  elsewhere (the sim was interrupted)
        
        #so you get back two iters one inside another
        #you could get back a certain line number an atom block very fast
        
    def idlinesbyatomid(self,atomblock):#useless
        #block has unique atoms
        #somehow the spc arg makes if faster
#        return (int(aline.split(' ')[self.attribs['id']]) \
#            for aline in atomblock)
    #MUCH faster!
        return (int(str.split(aline,' ')[self.attribs['id']]) \
            for aline in atomblock)
    
    def readnframes(self,n):
        """returns stacked frame info. breaks at diff no. of atoms"""
        fi=self.nextframeinfoi()
        tss=[];boxes=[];natoms=[];atomsblocks=[]
        #return [fi.next() for i in xrange(n)]
        def appender(af):
            tss.append(af[0]);boxes.append(af[1]);natoms.append(af[2]);
            atomsblocks.extend(af[3])
        
        ctr=0;lastna=None
        if self.broken!=None:ctr=1
        while ctr<n:
            if self.broken!=None:
                appender(self.broken)
                self.broken=None
            try:
                af=(fi.next())
                if lastna != None and af[2]!=lastna :
                    self.broken=af;break;
                appender(af)
                ctr+=1;lastna=af[2]
            except StopIteration:
                if ctr!=0:print ctr,'timesteps read';break
                else: raise StopIteration
        return {'timesteps':iter(tss),'boxbounds':iter(boxes),'natoms':iter(natoms)
        ,'atomsblocks':iter(atomsblocks)}

#    def aggregate(self,frames):
#        tss=(fi[0] for fi in frames)
#        boxes=(fi[1] for fi in frames)
#        natoms=(fi[2] for fi in frames)
#        atomblocks=([fi[3] for fi in frames])
#        stackedatomblocks=(item for sublist in atomblocks for item in sublist)
#        return tss,boxes,natoms,stackedatomblocks
               
#        tss=[];boxes=[];natoms=[]
#        atomblocks=[]
#        for aframe in frames:
#            tss.append(aframe[0]);boxes.append(aframe[1]);natoms.append(aframe[2])
#            atomblocks.extend(aframe[3])
#        return {'timesteps':iter(tss),'boxbounds':iter(boxes),'natoms':iter(natoms)
#        ,'atomblocks':iter(list(atomblocks))}#invoke the inside by list iter 1st

#todo interrupt seq if no. of atms changes. throw exception        
    def binatomlines(self,atomsblock):#in order
        """in:sequence of lines of atoms. out: dict of atoms's lines"""
        #apply numpy sort?
         #CRITICAL LOOPS here.
         #str.split (instead of aline.split) and ' ' and spec ' ' makes it faster
        sls=[str.split(aline,' ') for aline in atomsblock] #last item is \n
        #..don't want it
        lineids=[int(aline[self.attribs['id']]) for aline in sls]
        #nice
#        uids=frozenset(lineids) #could be skipped
#        #lesson: if lineids in a gen. it's consumed by frznst so i can only use
#        #it once but i need to use it more than once
#        store=dict(zip((uids),[[] for al in xrange(len((uids)))]))
#        #[[]]*len(uids))) #doesn't work!
#        def putinstore(aid,al):store[aid].append(al);return True#0 cheaper than None?
#        #map(putinstore,(lineids),(atomsblock))
#        map(putinstore,(lineids),sl)
        #better
        store=collections.defaultdict(list)
        for aid,asl in itertools.izip(lineids,sls): store[aid].append(asl)
        return store

    def extractvectorsfromlines(self,splitatomlines):
        """input one atom's lines"""
        vecsi=self.attribs.values();vecsi.remove(self.attribs['id'])
        sls=[operator.itemgetter(*vecsi)(al) for al in splitatomlines]
        #return numpy.transpose(sls) #makes it slow
        return dict(zip(self.attribs.keys(),zip(*sls))) #much faster! why? idk.
        #return zip(*sls)
        #return map(None,*sls) #slow
    
    def extractnext(self,n):
        framesdata=self.readnframes(n)
        #anatom,itsvecs
        binals=self.binatomlines(framesdata.pop('atomsblocks')) #should be no hit in pop
        eachatomwitsvectors=( (anatom, self.extractvectorsfromlines(binals[anatom])) \
            for anatom in binals.keys() )
        framesdata.update({'vecs':[eachatomwitsvectors]})
        return framesdata
        #want to map this op
#    def extractnext(self,atomid,colname):#vector index
#        vi
#        for ats in self.nexttimestep:
#        sl=str.split(self.scantonextatomline(atomid))
#        return sl[vi]
        #return operator.itemgetter(*vi)(sl)
    
    def user_extractallvectors(self,chunk=1000):
        while True:
            fd=self.extractnext(chunk)
            yield fd
        




#index is 3 make :3 len of atomidstr+spc
#[float(aline.split()[3]) for aline in d2vt.dumpfile if '10 ' in aline[:3]]
 #%prun d2vt.binatomlines(d2vt.aggregate([d2vt.nextframeinfoi().next()]*100)['atomblocks'])
#%prun list(d2vt.extractallvectors(400))