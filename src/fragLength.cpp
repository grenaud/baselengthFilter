#include <iostream>
#include <vector>
#include <cstring>
#include <cstring>
#include <set>
#include <ctype.h>
#include <stdlib.h>
#include <sys/mman.h>


#include <sys/types.h>

#include <sys/stat.h>
#include <fcntl.h>
#include "utils.h"

extern "C" {
    //#include "tabix.h"
    //#include "bam.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "bam.h"

#include "samtools.h"
#include "sam_opts.h"
#include "bedidx.h"

}


#define bam_is_reverse(b)     (((b)->core.flag&BAM_FREVERSE)    != 0)
#define bam_is_unmapped(b)    (((b)->core.flag&BAM_FUNMAP)      != 0)
#define bam_is_paired(b)      (((b)->core.flag&BAM_FPAIRED)     != 0)
#define bam_is_propaired(b)   (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam_is_read1(b)       (((b)->core.flag&BAM_FREAD1)      != 0)

#define bam_is_qcfailed(b)    (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_rmdup(b)       (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_sec(b)         (((b)->core.flag&BAM_FSECONDARY)        != 0)
#define bam_is_supp(b)        (((b)->core.flag&BAM_FSUPPLEMENTARY)    != 0)

#define bam_is_failed(b)      ( bam_is_qcfailed(b) || bam_is_rmdup(b) || bam_is_sec(b) || bam_is_supp(b) )

#define bam_mqual(b)          ((b)->core.qual)
#define bam_isize(b)          ((b)->core.isize)
#define bam_lqseq(b)          ((b)->core.l_qseq)




using namespace std;



// inline void minLFiltercout(const int32_t & l,const int32_t & m){
//     if( l>=m ){ 
// 	cout<<l<<endl; 
//     } 
// }

int main (int argc, char *argv[]) {

    bool onlyMapped    =false;
    bool onlyPP        =false;
    bool uncompressed  =false;
    bool filterQCfailed=false;
    
    int32_t  minLength =35;
    
    string usage=string(""+string(argv[0])+" <options>  [in BAM file]"+
			"\nThis program reads a BAM file and produces another with reads with a certain length cutoff\n"+
			//			"\nThis program reads a BAM file and produces another with reads with a certain length cutoff\n"+
			"\n"+

			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+

			"\n\n\tOther options:\n"+
			"\t\t"+"-l\t\t\tMinimum length of fragment to produce   (Default: "+booleanAsString( minLength )+")\n"+
			"\t\t"+"-m\t\t\tRequire the reads to be mapped          (Default: "+booleanAsString( onlyMapped )+")\n"+
			"\t\t"+"-p\t\t\tRequire the reads to be properly paired (Default: "+booleanAsString( onlyPP )+")\n"+
			"\t\t"+"-u\t\t\tProduce uncompressed BAM output         (Default: "+booleanAsString( uncompressed )+")\n"+
			"\t\t"+"-q\t\t\tFilter QC failed reads                  (Default: "+booleanAsString( filterQCfailed )+")\n"+
			
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    
    for(int i=1;i<(argc-1);i++){ //all but the last 3 args

        if(string(argv[i]) == "-l"  ){
            minLength=destringify<int32_t>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-m"  ){
            onlyMapped=true;
            continue;
        }

        if(string(argv[i]) == "-u"  ){
            uncompressed=true;
            continue;
        }

        if(string(argv[i]) == "-q"  ){
            filterQCfailed=true;
            continue;
        }

        if(string(argv[i]) == "-p"  ){
            onlyPP=true;
            continue;
        }

	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    
    string bamfiletopen    = string( argv[ argc-1 ] );
    string bamfiletwrte    = "/dev/stdout";
	
    
    samFile  *fpBAMin;
    samFile  *fpBAMout;
    bam1_t    *b;
    bam_hdr_t *h;
    
    fpBAMin  = sam_open_format(bamfiletopen.c_str(), "r", NULL);
    if(uncompressed)
	fpBAMout = sam_open_format(bamfiletwrte.c_str(), "wb0", NULL);
    else
	fpBAMout = sam_open_format(bamfiletwrte.c_str(), "wb", NULL);
    
    if(fpBAMin == NULL){
	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
	return 1;
    }
    
    if(fpBAMout == NULL){
	cerr << "Could not open output BAM file"<< bamfiletwrte << endl;
	return 1;
    }

    h = sam_hdr_read(fpBAMin);
    if(h == NULL){
	cerr<<"Could not read header for "<<bamfiletopen<<endl;
	return 1;
    }

    if(sam_hdr_write(fpBAMout, h) < 0){
	cerr<<"Could not write header for "<<bamfiletwrte<<endl;
	return 1;	
    }


    
    b = bam_init1();
    while(sam_read1(fpBAMin, h, b) >= 0){
	if(filterQCfailed)
	    if(bam_is_failed(b) )
		continue;
	
	if(bam_is_unmapped(b) ){
	    //if the read is unmapped and we only consider mapped reads
	    if(onlyMapped){  continue; }
	}else{//read is mapped	    
	    //we accept
	}

	//if(b->core.l_qseq < minLength) continue;
	bool ispaired    = bam_is_paired(b);
	//bool isfirstpair = bam_is_read1(b);
	
	if(ispaired){	    
	    //    if( isfirstpair   ){
	    if( !bam_is_propaired(b) ){
		if(onlyPP){  continue; }//if not properly paired, skip		    
	    }
	    
	}else{
	    //minLFiltercout(bam_lqseq(b),minLength);
	}

	if(bam_lqseq(b) < minLength) continue;
	
	if (sam_write1(fpBAMout, h, b) < 0) {    cerr<<"ERROR: Cannot write record to "<<bamfiletwrte<<endl; exit(1);}

    }
    
    bam_destroy1(b);
    
    if(sam_close(fpBAMin)<0){  cerr<<"Cannot close file: "<<bamfiletopen<<endl; }
    if(sam_close(fpBAMout)<0){ cerr<<"Cannot close file: "<<bamfiletwrte<<endl; }
   
    return 0;
}

