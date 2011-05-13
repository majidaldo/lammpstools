##Author: Majid S. al-Dosari
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

progdesc=\
"""
converts lammps dump files to hdf5 database. can resume if interrupted.

"""

import os
import sys

from dumps2vecs import dumps2vecs

import h5py
import numpy
class dumps2hdf5o(object):
    
    #filing is attrib/atom/vector from tuple description
    #timestep index
    def __init__(self,dumpreader,hdf5file):
        self.dbfilen=hdf5file
        self.dumpreader=dumpreader
        self.hdb=h5py.File(os.path.join(os.path.split(hdf5file)[0],'tmpdb')
        ,'a') #holding db
        #self.db=h5py.File(hdf5file,'a')
        self.attribdtype={}
        return
    def _heir2tuple(self,heir):
        return tuple([ab for ab in str.split(heir,'/') if ab is not ''])
    def _tuple2hier(self,tpl):
        if type(tpl)==str:return tpl
        path='/'
        for agroup in tpl: path=path+str(agroup)+'/'
        return path
    def returnrsrc(self,filing):
        return self.hdb[self._tuple2hier(filing)]
    def autodatatype(self,dataelem):
        """expecting strings from dump reader"""
        if type(dataelem)!=str: return type(dataelem)
        if '.' in dataelem: return 'f4'#float
        else: return 'i4'#int        
    def getattribdtype(self,attrib,dataelem):
        if attrib in self.attribdtype:return self.attribdtype[attrib]
        else:
            self.attribdtype.update({attrib:self.autodatatype(dataelem)})
            return self.attribdtype[attrib]
    
#    def createnewds(self,filing,arraydata):
#        ms=list(arraydata.shape)
#        ms[0]=None #only makes sense to append to 1st axis
#        self.hdb.create_dataset(filing,data=arraydata,maxshape=ms)
#        return self.hdb[filing]
    
    def appender1(self,arraywheader):#a dict w/ one item
        """input has to be array. {hdrforfiling (posix style): arraynotscalar}"""
        #filing is attrib/atom/vector 
        #check if index already there
        if len(arraywheader.values()[0])==0:return #ie nothing to add
        if type(arraywheader.keys()[0])!=str:#expect it to be a list
            filing=self._tuple2hier(arraywheader.keys()[0])
        else: filing =arraywheader.keys()[0] #posix style slashing
        arr=arraywheader.values()[0]
        if type(arr)!=numpy.array: arr=numpy.array( arr )
        if filing not in self.hdb:
            ms=list(arr.shape)
            ms[0]=None #only makes sense to append to 1st axis
            self.hdb.create_dataset(filing,data=arr,maxshape=ms
            ,chunks=True)#, compression='lzf') #can i create a blank?
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
                self.hdb.create_dataset(filing,data=pullout,maxshape=ms
                ,chunks=True)#, compression='lzf')
                self.hdb[filing].resize(arr.shape[0]+oldshape[0],axis=0)
                self.hdb[filing][oldshape[0]:]=arr
        return
    def appender(self,dictofarrayswheaders):
        for ahdr,anarray in dictofarrayswheaders.iteritems():
            self.appender1({ahdr:anarray})
        return        
    def appendifnotalreadyin(self,filing,tsarray):#for unordered but slow-er
        if self._tuple2heir(filing) in self.hdb:
            dbrsrc=self.hdb[self._tuple2hier(filing)]
            tsvecf=[apiece for apiece in tsarray if apiece[0] not in dbrsrc['ts']]
        else: tsvecf=tsarray #it's in there
        tsvecf=numpy.array(tsvecf #crap
            ,dtype=[('ts',int),('vals',self.getattribdtype(tsvecf[0][1]))])
        self.appender({filing:tsvecf })
        return
        
    def makeattribarray(self,tsi,attrib,sequence):
        #tsvec=zip(tsi,sequence) #zip was a bottleneck
        aa= numpy.zeros(len(tsi)
        ,dtype=[('ts',int),('vals',self.getattribdtype(attrib,sequence[0] ))])
        aa['ts']=tsi;aa['vals']=sequence
        return aa
        
        
    def resumefrom(self,under='/'):#for resuming
    #so then i can just scan to the ts in the dump
        tsis=self.returntsis(under=under)
        return min([max(atsi) for atsi in tsis])
    
    def atomsdumpi2hdfi(self,atomsdumpi,prefiling=[]):
        for stuff in atomsdumpi:
            tsi=list(stuff['timesteps']);na=list(stuff['natoms']);bb=list(stuff['boxbounds'])
            vecs=list(stuff['vecs'])
            vl=len(tsi)
            
            nar=numpy.zeros(vl, dtype=[('ts',int),('vals',int)])
            nar['ts']=tsi;nar['vals']=na
            self.appender1({tuple(prefiling+['natoms']):nar  } );del nar
#            self.appender1({tuple(prefiling+['timsteps']):numpy.array(zip(tsi,tsi)
#                ,dtype=[('ts',int),('vals',int)])    } )
            bbar=numpy.zeros(vl, dtype=[('ts',int),('vals','f4',(3,2))] )
            bbar['ts']=tsi;bbar['vals']=bb
            self.appender1({tuple(prefiling+['boxbounds']):bbar  });del bbar
            for anatomsstuff in vecs:
                aid=anatomsstuff[0];atomsdatadic=anatomsstuff[1]
                for anatomattrib,sequence in atomsdatadic.iteritems():                
                    filing=self._tuple2hier( tuple(prefiling+[anatomattrib]+[aid]))
                    #yield tsi,anatomattrib,sequence                    
                    tsvec=self.makeattribarray(tsi,anatomattrib,sequence)
                    self.appender1({filing:tsvec})
            yield tsi
    

    def sort(self, filing,by=['ts']):
        rsrc=self.returnrsrc(filing)
        #if type(rsrc)!=h5py.Dataset:raise Exception
        tosort=numpy.array(rsrc)
        tosort=numpy.sort(tosort,order=by)
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
        
    def uniquets(self,tsarray,tsname='ts'):
        utsi=numpy.unique(tsarray[tsname],return_index=True)[1]
        return tsarray[utsi]
    def removerepeatedts(self,adsn):#dataset name
        rrts=self.uniquets(self.hdb[adsn][:])
        # todo: don't do anything if the same?
        self.appender1({'/tmp/'+adsn:rrts})
        del self.hdb[adsn]
        fng=(self._heir2tuple(adsn))
#        self.hdb.copy('/tmp/toren',self._tuple2hier(fng[:-1])
#        ,name=self._tuple2hier(fng[-1]))#a sort of 'rename'
        self.appender1({fng:rrts})
        del self.hdb['/tmp/'+adsn]
        return rrts
        
    def removerepeatedtsindb(self,under='/',tsname='ts'):#useless
        under=self._tuple2hier(under)
        nms=[]
        def returndss(name,obj):
            if type(obj)!=h5py.Dataset \
            or 'tmp' in name:
                return None
            else: nms.append(name)
        
        g=self.returnrsrc(under)
        g.visititems(returndss)
        for adsn in nms:
            self.removerepeatedts(adsn)
        return
    
    def putinnewdb(self,under='/',tsname='ts'):
        """removes imperfections and compacts the db"""
        under=self._tuple2hier(under)
        nms=[]
        def returndss(name,obj):
            if type(obj)!=h5py.Dataset:
                return None
            else: nms.append(name)
        newdb=dumps2hdf5o(None,self.dbfilen)
        newdb.hdb.close()
        newdb.hdb=h5py.File(self.dbfilen,'a') #no tmp
        self.__init__(self.dumpreader,self.dbfilen)
        
        g=self.returnrsrc(under)
        g.visititems(returndss)
        
        for adsn in nms:
            norepeateddata=self.uniquets(self.hdb[adsn][:])
            if adsn not in newdb.hdb:
                newdb.appender1({adsn:norepeateddata['vals']})
            del self.hdb[adsn]
            #check this todo
            if nms[-1]==adsn:
                if 'timesteps' not in newdb.hdb:
                    newdb.appender1({'timesteps':norepeateddata[tsname]})
        return
    
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
    

        
pkwargs={'hdf5file': os.path.join(os.curdir,'atoms.hdf5')
,'chunk':1000
,'dumpdir':os.curdir,'dumpmatch':'*.dump'
,'deldumpfiles':False}
import glob
def dumps2hdf5(hdf5file=pkwargs['hdf5file']#output
,chunk=pkwargs['chunk'] #writing args
,dumpdir=pkwargs['dumpdir'],dumpmatch=pkwargs['dumpmatch']#dump reader args
,deldumpfiles=pkwargs['deldumpfiles']):
        
    dr=dumps2vecs(dumpdir=dumpdir,dumpmatch=dumpmatch)
    d2h=dumps2hdf5o( dr ,  hdf5file  )

    #if new db hasn't been created, means it's still reading/writing
    if not os.path.exists(d2h.dbfilen):   
        try:#route if resume needed
            rf=d2h.resumefrom() #this fails if there is no db
            #so it means that reading should continue
            while rf!=d2h.dumpreader.nexttimestep():continue
            print 'resuming reading from timestep', rf
        except:pass
        wg=d2h.atomsdumpi2hdfi(d2h.dumpreader.user_extractallvectors(chunk=chunk))
        print ''        
        for awi in wg:
            sys.stdout.write( "Completed writing upto timestep "
            +str(awi[-1])+'\b'*(len(str(awi[-1]))+32) )
            sys.stdout.flush()
        print ''
    
    #write to the new db
    print 'writing to final database'
    d2h.putinnewdb() #doesn't do anything if nothing there
    #del tmpdb
    df=d2h.hdb.filename
    d2h.hdb.close()
    try:print 'deleting temp database';os.remove(df)#)
    except: print 'could not delete temp db file'
    del df
    
    if deldumpfiles is True:
        print 'deleting dump files'
        od=os.path.abspath(os.path.curdir)
        os.chdir(dumpdir)
        dfs=glob.glob(dumpmatch)
        os.chdir(od);del od
        for adf in dfs: os.remove(adf)
        
    return h5py.File(hdf5file)


def cmdlinerun():
    from optparse import OptionParser

    parser = OptionParser(usage=progdesc)
    parser.add_option('--chunk',dest='chunk',type='int',default=pkwargs['chunk']
    ,help='how many timesteps to read at once. Default='+str(pkwargs['chunk'])+' seems best')
    parser.add_option('--dumpdir',dest='dumpdir',type='str',default=pkwargs['dumpdir']
    ,help='location of dump files. Default='+str(pkwargs['dumpdir'])   )
    parser.add_option('--dumpmatch',dest='dumpmatch',type='str',default=pkwargs['dumpmatch']
    ,help='dump file name pattern match. Default='+str(pkwargs['dumpmatch']))
    parser.add_option('--deldumpfiles',dest='deldumpfiles',action='store_true',default=pkwargs['deldumpfiles']
    ,help='delete dumpfiles after putting into database. Default='+str(pkwargs['deldumpfiles']))
    parser.add_option('--hdf5file',dest='hdf5file',type='str',default=pkwargs['hdf5file']
    ,help='destination database name. Default='+str(pkwargs['hdf5file']))
                      
    options, clargs = parser.parse_args()
    
    print options.dumpmatch
        
    dumps2hdf5(hdf5file=options.hdf5file
    ,chunk=options.chunk
    ,dumpdir=options.dumpdir,dumpmatch=options.dumpmatch
    ,deldumpfiles=options.deldumpfiles)
        
    return


if '__main__'==__name__:
    cmdlinerun()

