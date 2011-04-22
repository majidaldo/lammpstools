"""
use hd5py
subproc the 'main' so it collects results
"""

import os

olddir1=os.getcwd()
os.chdir(olddir1)
try:os.chdir('/home/aldosams/progs/lammps/tools/python/pizza')
except:os.chdir("\\\\129.59.197.166\\aldosams\\progs\\lammps\\tools\\python\\pizza")
from dump import dump as lammpsdump #instead of execfile("dump.py")
del olddir1


import fnmatch
import time



class dumpfilesreader(object):
    
    def __init__(self,dumpdir, colslist='all',atoms='all',matchpattern='*.dump'
        ,readbatch=10,stopson={}):        
        """better to have files to be fmt timestep.whatever"""
        self.stopson={'nonewfiles':False
        ,'timestep':False #ldump.time() #does it
        ,'numoffiles':False #num of read files, typically num of timesteps
        ,'extsignalflag':False} #a flag
        self.stopson.update(stopson)
        self.dumpdir=dumpdir;self.matchpattern=matchpattern
        self.readbatch=readbatch
        self.redfiles=[];self.redframes=[];self.okdumpfiles=[]
        self.waitfornewfiles()
        self.po=lammpsdump(os.path.join(dumpdir,matchpattern),0)
        if colslist=='all':
            self.po.next()
            self.cols=self.po.names.keys()
            self.cols.remove('id')#why would i want id?!
            #reinit
            self.po=lammpsdump(os.path.join(dumpdir,matchpattern),0)
        else: self.cols=colslist
        return
    
#    def readbatch(self):
#        """reads 10 files or timesteps whichever first"""
#        red=0
        
    
    def reinitdump(self,files):#todo w each new 'scan'?
        self.po.lammpsdump(files)
        
    def chkheadingsexist(self,adf):#a dump file
        #the headings should exist at least when before being read
        fo=open(adf)
        foundtimestep=False;foundnumatms=False;foundboxbounds=False;foundatoms=False
        #cont=True
        lcount=0 
        while foundtimestep==False\
        or foundnumatms==False\
        or foundboxbounds==False\
        or foundatoms==False:
            al=fo.readline()
            lcount+=1
            if 'TIMESTEP' in al:foundtimestep=True;continue
            if 'NUMBER' in al: foundnumatms=True;continue
            if 'BOX' in al: foundboxbounds=True;continue
            if 'ATOMS' in al: foundatoms=True;continue
            #but if
            if lcount>10:#shouldn't need to go down alot.
                fo.close()
                return False
        fo.close()
        return True
    def areallatomsthere(self,adf):#chk getatomis exception
        try:self.getatomids1(adf);return True
        except: return False
    def chkdumpfile1(self,adf):
        if self.chkheadingsexist(adf)==True:
            if self.areallatomsthere(adf)==True: return True
            else: return False
        else: return False
    def chkdumpfile(self,adf):
        if adf in self.okdumpfiles: return True
        else:
            dfc=self.chkdumpfile1(adf)
            if dfc==True: self.okdumpfiles.append(adf);return dfc
            else: return dfc #should be false
        
    def waitfornewfiles(self):
        while self.aretherenewfiles()==False: time.sleep(3)
        return

    def returnnewfiles(self):
        #so even if file was deleted it won't be picked up again
        l= list(frozenset(self.getdumpfiles())-frozenset(self.redfiles))
        l.sort()        
        l=self.chkdumpfiles(l[:self.readbatch])
        l.sort()
        return l
        
    def aretherenewfiles(self):
        if len(self.returnnewfiles())>0: return True
        else: return False
        
        #rename readbatch
    def readfilesbatch(self):#assuming there will be, doesn't stop
        nfs=self.returnnewfiles()[:self.readbatch]
        #clip extra files in a batch
        if self.stopson['numoffiles']!=False\
        and len(nfs)+len(self.redfiles)>self.stopson['numoffiles']:
            nxtra=len(nfs)+len(self.redfiles)-self.stopson['numoffiles']
            nfs=nfs[:-nxtra]
        self.po=lammpsdump(self._list2strwspcs(nfs))
#       self.po=lammpsdump(self._list2strwspcs(nfs),0)
#        red=0
#        while red<self.readbatch:#should be able to read a huge file
#        #while it's being written
#            if self.po.next()==-1: break
#            else: red+=1
        if self.stopson['timestep']!=False:
            self.po.tselect.test("$t <= "+str(self.stopson['timestep']))
            self.po.delete() #erases extras in case of last iteration
        urf=self.returnunreadfiles(nfs)
        self.redfiles.extend(list(frozenset(nfs)-frozenset(urf)))
        return #stderr if no dump files
    
    def returnunreadfiles(self,dfs):#dumpfiles
        """to take care of cases when a frame is being written while being read"""
        readframes=self.po.time()
        readfiles=[]
        for adf in dfs:
            for af in readframes:
                if str(af) in adf: readfiles.append(adf)
        return list(frozenset(dfs)-frozenset(readfiles))
    
    def evalstop(self):
        stopson=self.stopson
        if stopson['extsignalflag']==True:return True
        if stopson['nonewfiles']!=False:
            if self.aretherenewfiles()==False:return True
        if stopson['timestep']!=False:
            for afn in self.redfiles:
                if str(stopson['timestep']) in afn\
                and self.aretherenewfiles()==False: return True
        if stopson['numoffiles']!=False:
            if len(self.redfiles)>=stopson['numoffiles']:return True
        #else
        return False
    
    def _list2strwspcs(self,lst):
        lst=str(lst)
        lst=lst.replace('[','');lst=lst.replace(']','')
        lst= lst.replace("'",'')
        lst= lst.replace("\\\\",'\\') #why do i have to do this?    
        return lst.replace(',','')
    
    #make a dump file iterator w/ glob litmit it w/ readbatch
    def getdumpfiles(self,matchpattern='*.dump'):#kwargs update thing
        dumpdir=self.dumpdir;matchpattern=self.matchpattern
        walker=os.walk(dumpdir)
        files=walker.next()[2];del walker
        fh=fnmatch.filter(files,matchpattern);del files
        fullfns=[os.path.join(dumpdir,afile) for afile in fh]
        #fullfns=fullfns[:self.readbatch]
        return fullfns
        
    def chkdumpfiles(self,dfs):
        #screen for completeness
        icfs=[af for af in dfs if False==self.chkdumpfile(af)]
#        for af in fullfns:
#            if True==self.chkdumpfile(af): icfs.append(af)
        return list(frozenset(dfs)-frozenset(icfs))
    
    def getcolnames(self):#,put in a/1 dumpfile arg?
        #useless? pizza.names after a pizza.next
        """assume group in the same dir has same cols. the dump file should be
        small"""
        anydf=self.getdumpfiles()[0]#take 1st one i get
        fo=open(anydf)
        cont=True
        while cont==True:
            al=fo.readline()
            if 'ITEM: ATOMS' in al: #ie line is: ITEM: ATOMS id type vx vy vz 
                cols=al.split()[3:]
                cont=False
            else: pass
        fo.close()
        return cols
    
    def getatomids(self):
        anydf=self.getdumpfiles()[0]#take 1st one i get
        return self.getatomids1(anydf)
    def getatomids1(self,adumpfile):
        #can be used to chk completeness of file
        anydf=adumpfile
        fo=open(anydf)
        cont=True
        atomids=[]
        #get no. of atoms
        while cont==True:
            al=fo.readline()
            if 'ITEM: NUMBER OF ATOMS' in al:
                #read the no. in the next line
                natoms=int(fo.readline())              
                cont=False
        if natoms==0: return atomids
        #then read to titles
        cont=True
        while cont==True:
            al=fo.readline()
            if 'ITEM: ATOMS' in al: cont=False
        #so the next line should have an atom
        cont=True;atomcount=0
        while atomcount<natoms:
            try:al=fo.readline()
            except EOFError: raise Exception, '1st frame incomplete'
            try:anatomid=al.split()[0]
            except: raise Exception, '1st frame incomplete'
            atomcount+=1;atomids.append(int(anatomid))
        fo.close()
        return atomids
        
    def returnatomvecs(self,atmid,veclist):#atmname, vectorname
    #used after a read
        #vectors=self.po.fasteratom(atmid,*veclist)
        vectors=self.po.atom(atmid,*veclist)
        if len(veclist)==1:vectors=[vectors]
        vectors=dict(zip(veclist,vectors))
        return self.po.time(),vectors #order, vecs
    def returnatomsvecs(self,atomids,veclist):
        o={}   
        for anatom in atomids:
            ao=self.returnatomvecs(anatom, veclist)
            for anattrib,avec in ao[1].iteritems():
                o.update({(anattrib,anatom):avec})
        return ao[0],o
        
    def dumpreaderi(self):
        while self.evalstop()==False:
            self.waitfornewfiles() 
            #if len(self.redfiles)>1: self.waitfornewfiles() #that it,
            #..the intention is to write one frame to a file instead of one
            #huge file
            yield self.readfilesbatch()
            
    def applyf(self,afunc,itsargs):
        di=self.dumpreaderi()
        for ani in di:
            yield afunc(*itsargs)
            
    def atomsdumpi(self,atomslist,attribslist):
        d2hdfi=self.applyf(self.returnatomsvecs,(atomslist,attribslist))
        return d2hdfi
            
    
        
import h5py
import scipy
class dump2hdf5(object):
    
    #filing is attrib/atom/vector from tuple description
    #timestep index
    def __init__(self,dumpreader,hdf5file):
        self.dumpreader=dumpreader
        self.hdb=h5py.File(hdf5file,'a')
        return
    
    def _tuple2hier(self,tpl):
        if type(tpl)==str:return tpl
        path='/'
        for agroup in tpl: path=path+str(agroup)+'/'
        return path
    def returnrsrc(self,filing):
        return self.hdb[self._tuple2hier(filing)]
    
#    def createnewds(self,filing,arraydata):
#        ms=list(arraydata.shape)
#        ms[0]=None #only makes sense to append to 1st axis
#        self.hdb.create_dataset(filing,data=arraydata,maxshape=ms)
#        return self.hdb[filing]
    
    def appender1(self,arraywheader):#a dict w/ one item
        """input has to be array. {hdrforfiling: arraynotscalar}"""
        #filing is attrib/atom/vector 
        #check if index already there
        if len(arraywheader.values()[0])==0:return #ie nothing to add
        if type(arraywheader.keys()[0])!=str:#expect it to be a list
            filing=self._tuple2hier(arraywheader.keys()[0])
        else: filing =arraywheader.keys()[0] #posix style slashing
        arr=scipy.array(arraywheader.values()[0])
        if filing not in self.hdb:
            ms=list(arr.shape)
            ms[0]=None #only makes sense to append to 1st axis
            self.hdb.create_dataset(filing,data=arr,maxshape=ms)#, compression='lzf') #can i create a blank?
        else:
            oldshape=self.hdb[filing].shape #a copy
            try:
                self.hdb[filing].resize(arr.shape[0]+oldshape[0],axis=0)
                self.hdb[filing][oldshape[0]:]=arr
            except TypeError:#if the data was created w/o maxshape=None
                pullout=self.hdb[filing][:]
                del self.hdb[filing]
                ms=list(arr.shape)
                ms[0]=None #only makes sense to append to 1st axis
                self.hdb.create_dataset(filing,data=pullout,maxshape=ms)#, compression='lzf')
                self.hdb[filing].resize(arr.shape[0]+oldshape[0],axis=0)
                self.hdb[filing][oldshape[0]:]=arr
        return
    def appender(self,dictofarrayswheaders):
        for ahdr,anarray in dictofarrayswheaders.iteritems():
            self.appender1({ahdr:anarray})
        return
        
    def atomsdumpi2hdfi(self,atomsdumpi):
        for stuff in atomsdumpi:
            #[0]is timestep index,   [1] is vecs          
            #todo here should go chking for existing ts,aid,attrib
            #self.appender({'tsi':stuff[0]})
            #stuff2={} #to put vectors with (ts,value)
            for afiling,avec in stuff[1].iteritems():
                #stuff2.update({afiling:zip(stuff[0],avec)})
                tsvec=zip(stuff[0],avec)
                #chk if it's not in already
#                tsvec2=tsvec[:]
#                for apiece in tsvec:
#                    if self.checkifexists1(afiling,apiece)==True:
#                        tsvec2.remove(apiece)
                if self._tuple2hier(afiling) in self.hdb:
                    #print tsvec[0]
                    #print tsvec[0] in self.hdb[self._tuple2heir(afiling)]
                    #^^^doesn't work?!
                    dbr=self.hdb[self._tuple2hier(afiling)] #db rsrc
                    tsvecf=[apiece for apiece in tsvec \
                        if apiece[0] not in dbr['ts']]
                else: tsvecf=tsvec
                #filter(self.checkifexists1,tsvec)
                tsvecf=scipy.array(tsvecf,dtype=[('ts',int),('vals','f4')])
                self.appender({afiling:tsvecf })
            #self.appender(stuff[1])
                yield afiling ,tsvecf #what went in

    def atomsvecs2dhfi(self,atomslist,vecslist):#todo
        di=self.dumpreader.atomsdumpi(atomslist,vecslist)
        #then put it in atomsdumpi2hdfi
        return self.atomsdumpi2hdfi(di)
    
#    def checkifexists1(self,filing, datapiece):
#        if self._tuple2heir(filing) and datapiece in self.hdb[self._tuple2heir(filing)]:
#            return True
#        else: return False
    def sort(self, filing,by=['ts']):
        rsrc=self.returnrsrc(filing)
        #if type(rsrc)!=h5py.Dataset:raise Exception
        tosort=scipy.array(rsrc)
        tosort=scipy.sort(tosort,order=by)
        self.hdb[self._tuple2hier(filing)][:]=tosort
        return tosort #ValueError if unknown field    
        
    def sortall(self,under='/',by=['ts']):
        under=self._tuple2hier(under)
        def f(name,obj):
            if type(obj)!=h5py.Dataset:return None#raise Exception
            try:self.sort(name,by=by)
            except: pass
            return
        g=self.returnrsrc(under)
        v=g.visititems(f)
        return v
    
    def returntsis(self,under='/',tsname='ts'):
        under=self._tuple2hier(under)
        tsis={}
        def f(name,obj):
            if type(obj)!=h5py.Dataset:return None#raise Exception
            try:
                tsi=obj[tsname] #to see if it has an index
                try:tsis[frozenset(tsi)].append(name)#see if it's there already
                except:tsis.update({frozenset(tsi):[name]})
            except: pass
            return
        g=self.returnrsrc(under)
        g.visititems(f)
        return tsis
        #if just one item in this then OK

    #make it begin from the shortest tsi. todo so need to add
    # a beginning this for dump reader
        

        #todo: del files
def dump2hd5program(dbpath,dumpdir
,atomids='all',attribs='all'
,readbatch=10
,stopson={'nonewfiles':True}):
    d2vs=dumpfilesreader(dumpdir,stopson=stopson,readbatch=readbatch)
    d2db=dump2hdf5(d2vs,dbpath)
    #progress=[]
    
    if atomids=='all':atomids=d2vs.getatomids()
    if attribs=='all':attribs=d2vs.getcolnames()
    #just to create a file if it doesn't exist
    #if os.path.isfile(dbpath):fh=open(dbpath,'w');fh.close()
    it=d2db.atomsvecs2dhfi(atomids,attribs)
    #return it
    for ani in it:
        #if ani[0] not in progress: progress.append(ani[0]); print ani[0]
        #print d2vs.evalstop()
        continue
    print 'sorting all vectors...'
    d2db.sortall()
    
    tsis=d2db.returntsis()
    if len(tsis)!=1: return tsis #you have a problem
    return d2db.hdb
    #select all atoms to dump or partial
    # creat file if doesn't exist
    #saves state, can resume
    #ctrl merging diff passes
    #deletes files new idea ..if redfile open, delete
    # optionally for cases where each file is a 
    # make a case for when each file is a frame. write method
    #close db file

#another proc deletes
#will need to sort in the db maybe
    #then matrix atmids and vns


testdir="\\\\129.59.197.166\\aldosams\\anontplab\\research\\yag\\runtypes\\eqvib\\0\\dump"
bdf="\\\\129.59.197.166\\aldosams\\anontplab\\research\\yag\\runtypes\\eqvib\\0\\"
dbf="\\\\129.59.197.166\\aldosams\\anontplab\\research\\yag\\runtypes\\eqvib\\0\\d.hdf5"

d2vt=dumpfilesreader(testdir
,stopson={'nonewfiles':True}
##,stopson={'numoffiles':15})
##,stopson={'nonewfiles':True})
##{'numoffiles':11})
##,stopson={'timestep':5000}
)
#d2ht=dump2hdf5(d2vt, dbf)

#just put it in hdfs5
#something that takes frames and stacks it into vecs
    
        
        #        fns=self.batchobj.user_getoutputvsinput(
#        params,listofpatterns=['*.dump']).values()[0]
        
    
    
    #todo supress stdout
        
    #todo dump reader its own class
    #eqviba.readatomattrib(eqviba.params[0],1,['vx'])
#    def readatomattrib(self,params,atomid,attribslist):
#        fns=self.batchobj.user_getoutputvsinput(
#        params,listofpatterns=['*.dump']).values()[0]
#        print fns
#        fns=str.strip(str(fns),'[]')
#        fns=fns.replace(',',' ');fns=fns.replace("'",'')
#        fns=fns.replace('\\\\','\\') #weird!
#        #return fns #does this work on linux?        
#        do=lammpsdump(fns,0)
#        av={}
#        curts=do.next()
#        while curts!=-1 and curts<100100: #for testing
#            print curts
#            aval=do.atom(atomid,*attribslist)
#            av.update({curts:aval})
#            do.tselect.one(curts);do.delete() #erase previous
#            curts=do.next()
#        return av
