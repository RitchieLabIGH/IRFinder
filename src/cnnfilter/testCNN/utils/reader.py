import json
import gzip

class ImageArray:
    def __init__(self, raw_line):
        raw_line=raw_line.split("\t")
        raw_name=raw_line[0].split(":")
        self.name=raw_name[0]+":"+raw_name[1]
        pos=raw_name[1].split("-")
        self.intron_name="{}:{}-{}".format(raw_name[0],int(pos[0])+15,int(pos[1])-15)
        self.strand=raw_name[2]
        self.region=json.loads(raw_line[1])
        self.is_valid= len(self.region) > 1
        if self.is_valid:
            self.is_valid=max([sum(i) for i in self.region]) > 0
    

class ImageArchive:
    
    def _open_file(self, fname):
        if fname.endswith(".gz"):
            return gzip.open(fname, "rt")
        else:
            return open(fname, "rt")
    
    def _count_lines(self, fname):
        tmp = self._open_file(fname)
        count=0
        for _ in tmp:
            count+=1
        tmp.close()
        return count
        
        
    def __init__(self, bed_file, array_file):
        self._len=self._count_lines(array_file)
        if bed_file != None:
            self._bed = self._open_file(bed_file)
            self._has_bed=True
        else:
            self._has_bed=False
        self._array= self._open_file(array_file)
        self._index=-1
        
            
    def __iter__(self):
        return self
    
    def __next__(self):
        self._index+=1
        if self._index < self._len:
            if self._has_bed:
                img_array=ImageArray(self._array.readline())
                bed_line=self._bed.readline().strip().split("\t")
                while bed_line[0] + ":" + bed_line[1] + "-" + bed_line[2] != img_array.intron_name :
                        bed_line=self._bed.readline().strip().split("\t")
                return bed_line , img_array
            else:
                return ["NA"], ImageArray(self._array.readline())
        else:
            raise StopIteration
    
    def __len__(self):
        return self._len
    
    def __del__(self):
        self.close()
        
    def close(self):
        self._array.close()
        if  self._has_bed:
            self._bed.close()
        
    def getIndex(self):
        return self._index
    
