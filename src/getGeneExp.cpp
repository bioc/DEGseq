#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

using namespace std;
class Exon;
int getCounts(string refFlatFile, vector<string> MapResultFiles, string outputFile, string Format, bool needSameStrand);
int getCountsForOneSample(vector<Exon> exons, map<string,int> maxExonLen, vector<string> geneNames,
                          string MapResultFile, string Format, bool needSameStrand,
                          map<string,pair<int,int> > & genesExp, map<string, int> & readsCount);
int getAllFiles(vector<string> FilesAndDirs, vector<string> & AllFiles);
int getExons(string refFlatFile, vector<Exon> &exons, vector<string> &geneNames, map<string,int>  & maxExonLen, map<string, int> & geneAveLen);
int getExonAnnotationFile2(char ** refFlatFile_str, char ** exonAnnotationFile);
static int rows = 0;
static double overlap_percent = 1;
static bool addtionCol = false;
static string path_sep = "\\";
extern "C" {
int getGeneExp(char ** refFlatFile_str, char ** MapResultFile_str, int * fileCount, char ** outputFile_str,
                 char ** Format_str, int * readlength_int, int * needSameStrand_int, double * overlapPercent)
{

  string refFlatFile = refFlatFile_str[0];
  vector<string> MapResultFiles;
  for(int i=0; i < fileCount[0]; i++ ){
      MapResultFiles.push_back(MapResultFile_str[i]);
  }
  vector<string> AllFiles;
  getAllFiles(MapResultFiles, AllFiles);
  string outputFile = outputFile_str[0];
  string Format = Format_str[0];
  int readlength = readlength_int[0];
  int needSameStrand = needSameStrand_int[0];
  overlap_percent = overlapPercent[0];
  if(Format==""){
     Format="bed";
  }
  if((refFlatFile.c_str()[0]=='"')&&(refFlatFile.c_str()[0]=='"')){
     refFlatFile=refFlatFile.substr(1,refFlatFile.size()-2);
  }
  if((outputFile.c_str()[0]=='"')&&(outputFile.c_str()[0]=='"')){
     outputFile=outputFile.substr(1,outputFile.size()-2);
  }

  if(Format == "eland"){
      if(readlength==0){
         readlength=32;
      }else{
      }
      char tmp[10];
      sprintf(tmp,"%d",readlength);
      Format= Format+tmp;
  }
  clock_t start = clock();
  Rprintf("Count the number of reads mapped to each gene.\n");
  Rprintf("This will take several minutes.\nPlease wait ...\n");
  R_FlushConsole();
  if(getCounts(refFlatFile, AllFiles, outputFile, Format, needSameStrand) < 0){
    Rprintf("There is something wrong!Please check...\n");
    R_FlushConsole();
    return -1;
  }
  clock_t end = clock();
  Rprintf("total used %f seconds!\n",((double)end - start)/CLOCKS_PER_SEC);
  R_FlushConsole();
  rows=0;
  overlap_percent=1;
  addtionCol = false;
  return 0;
}
}
extern "C" {
int getExonAnnotationFile(char ** refFlatFile_str, char ** exonAnnotationFile)
{
    return getExonAnnotationFile2(refFlatFile_str, exonAnnotationFile);
}
}

int isDir(const char * name){
  struct stat buff;
  char name2[200];
  strcpy(name2, name);
  int len = strlen(name2);
  if((name2[len-1] == '\\')||(name2[len-1] == '/')){
      name2[len-1] = '\0';
  }
  if(stat(name2, &buff) < 0)
     return -1;
  return S_ISDIR(buff.st_mode);
}
string get_file_name(const string file_name) {
  size_t i = file_name.find_last_of("//");
  if (i == string::npos){
      i = file_name.find_last_of("/");
      if(i == string::npos){
         i = file_name.find_last_of("\\");
      }
  }
  return file_name.substr(i+1);
}

int getAllFiles(vector<string> FilesAndDirs, vector<string> & AllFiles){
  AllFiles.clear();
  for(vector<string>::iterator it=FilesAndDirs.begin();it!=FilesAndDirs.end();it++){
      string name=*it;
      if((name.c_str()[0]=='"')&&(name.c_str()[0]=='"')){
        name=name.substr(1,name.size()-2);
      }
      char filename[200];
      strcpy(filename, name.c_str());
      int len = strlen(filename);
      if((filename[len-1] == '\\')||(filename[len-1] == '/')){
          filename[len-1] = '\0';
      }
      string dir_str = filename;
      if(isDir(filename) < 0){
          Rprintf("%s does not exist!\n", it->c_str());
	        R_FlushConsole();
          //return -1;
          continue;
      }
      if(isDir(filename)){
         DIR *dir;
         if (!(dir = opendir(filename))){
         }
         struct dirent *ent;
         while ((ent = readdir(dir))) {
            if(ent->d_name[0]=='.') continue;
            AllFiles.push_back(dir_str+path_sep+string(ent->d_name));
         }
      }else{
         AllFiles.push_back(name);
      }
  }
  Rprintf("SampleFiles:\n");
  for(vector<string>::iterator it=AllFiles.begin();it!=AllFiles.end();it++){
      Rprintf("\t%s\n",it->c_str());
  }
  R_FlushConsole();
  return 0;
}

void string2tokens(const string& str, vector<string> &tokens ,
                   const char delimiters = '\t', bool skip_empty = true) {

        int lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
        lastPos--;
        int pos = str.find_first_of(delimiters, lastPos+1);
        tokens.clear();
        while(true)
        {
            if (pos == -1)
            {
                if(lastPos+1>=(int)str.size()){
                   return;
                }
                tokens.push_back(str.substr(lastPos+1, str.size()-lastPos-1));
                return;
            }
            if (pos>lastPos+1)
            {
                tokens.push_back(str.substr(lastPos+1, pos - lastPos-1));
            }

            lastPos = pos;
            pos=str.find_first_of(delimiters, lastPos+1);
        }
        return;
}


class Exon{
    public:
       string geneName;
       string chr;
       int start;
       int end;
       int strand;
       Exon(string geneName, string chr, int start,int end,int strand=0)
               :geneName(geneName),chr(chr),start(start),end(end),strand(strand){ }
       inline bool operator<(const Exon &e) const {
           if(chr != e.chr){
               return chr < e.chr;
           }
           if(start != e.start){
               return start < e.start;
           }
           if(end != e.end){
               return start < e.start;
           }
           if(strand != e.strand ){
               return strand < e.strand; //we assume "+" < "-"
           }
           if(geneName != e.geneName ){
               return geneName < e.geneName;
           }
           return false;
       }

       inline bool operator==(const Exon &e) const {
           if(chr != e.chr){
               return false;
           }
           if(start != e.start){
               return false;
           }
           if(end != e.end){
               return false;
           }
           if(strand != e.strand ){
               return false;
           }
           if(geneName != e.geneName ){
               return false;
           }
           return true;
       }
};

ostream& operator<<(ostream& s, const Exon& exon) {
  string strand("-");
  if(exon.strand==0){
     strand="+";
  }
  return s<<exon.geneName<<"\t"<<exon.chr<<"\t"<<exon.start<<"\t"<<exon.end<<"\t"<<strand<<endl;
}

class Read{
    public:
       string chr;
       int start;
       int end;
       int strand;  //0:+   1:-
       int count;
       Read():chr("chr0"),start(0),end(0),strand(0),count(0){}
       Read(string chr, int start,int end,int strand):chr(chr),start(start),end(end),strand(strand),count(1){ }
       Read(string chr, int start,int end,int strand,int count):chr(chr),start(start),end(end),strand(strand),count(count){ }
       inline bool operator<(const Read &r) const {
           if(chr != r.chr){
               return chr < r.chr;
           }
           if(start != r.start){
               return start < r.start;
           }
           if(end != r.end){
               return start < r.start;
           }
           if(strand != r.strand ){
               return strand < r.strand; //we assume "+" < "-"
           }
           return false;
       }

       inline bool operator==(const Read &r) const {
           if(chr != r.chr){
               return false;
           }
           if(start != r.start){
               return false;
           }
           if(end != r.end){
               return false;
           }
           if(strand != r.strand ){
               return false;
           }
           return true;
       }

       inline bool operator<(const Exon &e) const {
           if(chr != e.chr){
               return chr < e.chr;
           }
           if(end <= e.start){
               return true;
           }
           return false;
       }

       inline bool operator>(const Exon &e) const {
           if(chr != e.chr){
               return chr > e.chr;
           }
           if(start >= e.end){
               return true;
           }
           return false;
       }

       inline int overlap(const Exon &e) const {
           if(chr != e.chr){
              return 0;
           }
           if(end <= e.start){                         //r1 r2 e1 e2
               return 0;
           }
           if(start >= e.end){                         //e1 e2 r1 r2
               return 0;
           }
           if((start >= e.start)&&(end <= e.end)){     //e1 r1 r2 e2
               return end-start;
           }
           if((start >= e.start)&&(end > e.end)){      //e1 r1 e2 r2
               return e.end-start;
           }
           if((start <= e.start)&&(end >= e.end)){     //r1 e1 e2 r2
               return e.end-e.start;
           }
           if((start <= e.start)&&(end < e.end)){     //r1 e1 r2 e2
               return end-e.start;
           }
           Rprintf("There is something wrong!\n");
      	   R_FlushConsole();
           return -1;
       }
};

ostream& operator<<(ostream& s, const Read& read) {
  string strand("-");
  if(read.strand==0){
     strand="+";
  }
  return s<<read.chr<<"\t"<<read.start<<"\t"<<read.end<<"\t"<<strand<<endl;
}

int getExons(string refFlatFile, vector<Exon> &exons, vector<string> &geneNames, map<string,int>  & maxExonLen, map<string, int> & geneAveLen){
    exons.clear();
    ifstream in(refFlatFile.c_str());
    if(!in) {
      Rprintf("cannot open input file %s\n", refFlatFile.c_str());
      R_FlushConsole();
      return -1;
    }
    map<string, vector<int> > geneLen;
    
    int count=0;
    while (!in.eof()) {
         //printf("\r%d",count);
         count++;
         char buffer[10000 + 1];
         in.getline(buffer, 10000);
         if (buffer[strlen(buffer) - 1] == '\r'){
             //buffer[strlen(buffer) - 1] = '\0';
         }
         int isoformLen=0;
         vector<string> blocks;
         string tmp = buffer;
         if(tmp.size()<10){
             continue;
         }
         string2tokens(string(buffer),blocks);
         if(blocks.size()<5){
             continue;
         }
         vector<string> begins;
         vector<string> end;
         string geneName = blocks[0];
         string chr=blocks[2];
         int strand=0;
         if(blocks[3]=="+"){
            strand = 0;
         }else if(blocks[3]=="-"){
            strand = 1;
         }else{
            Rprintf("Wrong refFlat format!\n");
      	    R_FlushConsole();
            return -1;
         }
         vector<string> exonSrt;
         vector<string> exonEnd;

         string2tokens(blocks[9],exonSrt,',');
         string2tokens(blocks[10],exonEnd,',');
         int iter = 0;

         for(vector<string>::iterator it=exonSrt.begin();it!=exonSrt.end();it++){
             int start = atoi(it->c_str());
             int end = atoi(exonEnd[iter].c_str());
             Exon exon(geneName,chr,start,end,strand);
             exons.push_back(exon);
             if(maxExonLen.count(chr) == 0){
                maxExonLen[chr] = end-start;
             }else{
                 if(maxExonLen[chr] < (end-start)){
                    maxExonLen[chr] = end-start;
                 }
             }
             isoformLen += end-start;
             iter++;
         }
         geneNames.push_back(geneName);
         if(geneLen.count(geneName) ==0){
         	  vector <int> tmp;
         	  tmp.push_back(isoformLen);
         	  geneLen[geneName] = tmp;
         }else{
            geneLen[geneName].push_back(isoformLen);
         }
    }
    //printf("\nBegin the sort!!!%d\n",exons.size());
    sort(exons.begin(),exons.end());
    vector<Exon>::iterator ptr1 = unique(exons.begin(),exons.end());
    exons.erase(ptr1,exons.end());

    sort(geneNames.begin(),geneNames.end());
    vector<string>::iterator ptr2 = unique(geneNames.begin(),geneNames.end());
    geneNames.erase(ptr2, geneNames.end());
    Rprintf("total %d unique genes\n",(int)geneNames.size());
    R_FlushConsole();
    
    
    for(vector<string>::iterator it=geneNames.begin();it!=geneNames.end();it++){
    	  int total_len = 0;
    	  for(vector<int>::iterator it2=geneLen[*it].begin();it2!=geneLen[*it].end();it2++){
    	  	  total_len += *it2;
    	  }
    	  geneAveLen[*it] = int(total_len/geneLen[*it].size());
    }
    return 0;
}


int findOverLapGenes(vector<Exon> & exons, Read & one_read, map<string,pair<int,int> > & genesExp, map<string,int> & maxExonLen, bool needSameStrand){
    int head=0;
    int tail=exons.size()-1;
    int pos;
    while(1){
      pos=(head+tail)/2;
      if(pos==head){
         pos=head+1;
         break;
      }
      if((!(one_read < exons[pos]))&&(one_read < exons[pos+1])){
         break;
      }else if(one_read < exons[pos]){
         tail = pos-1;
         continue;
      }else if(!(one_read < exons[pos+1])){
         head = pos;
         continue;
      }else{
         Rprintf("bug!\n");
	       R_FlushConsole();
         return -1;
      }
    }

    int readLen = one_read.end-one_read.start;

    vector<string> geneNames;
    double overLen = overlap_percent*readLen;
    if(overLen <= 1){
      overLen=1;
    }else if(overLen > readLen){
      overLen=readLen;
    }
    while(pos>=0){
        if(exons[pos].chr < one_read.chr){
           break;
        }
        if(exons[pos].start + maxExonLen[one_read.chr] < one_read.start){
           break;
        }
	if((needSameStrand==true)&&(exons[pos].strand != one_read.strand)){
	  pos--;
	  continue;
	}
        if(one_read.overlap(exons[pos])+0.001 >= overLen){
           geneNames.push_back(exons[pos].geneName);
        }
        pos--;
    }

    sort(geneNames.begin(),geneNames.end());
    vector<string>::iterator ptr = unique(geneNames.begin(),geneNames.end());
    geneNames.erase(ptr, geneNames.end());
    if(geneNames.size()==1){
        genesExp[(*geneNames.begin())].first++;
        genesExp[(*geneNames.begin())].second++;
    }else{
        for(vector<string>::iterator it = geneNames.begin();it!=geneNames.end();it++){
            genesExp[(*it)].second++;
        }
    }
    return 0;
}


void printResult(string outputFile, map<string,pair<int,int> > & genesExp, vector<string> geneNames){
    ofstream out(outputFile.c_str());
    if(!out) {
      cout<<string("cannot open output file ") + outputFile<<endl;
      exit(0);
    }
    out<<"\"geneName\"\t"<<"\"reads uniquely mapped to this gene\"\t"
            <<"\"reads mapped to this gene\"\t"<<"\"reads uniquely mapped to genome\""<<endl;
    for(vector<string>::iterator it=geneNames.begin();it!=geneNames.end();it++){
       out<<*it<<"\t"<<genesExp[*it].first<<"\t"<<genesExp[*it].second<<"\t"<<rows<<endl;
    }
}
void printResult2(string outputFile, vector<string> MapResultFiles, map<string, map<string, pair<int,int> > > & genesExp,
                  vector<string> geneNames, map<string, int> & readsCount){
    ofstream out(outputFile.c_str());
    if(!out) {
      Rprintf("cannot open output file %s \n", outputFile.c_str());
      return;
    }
    out<<"\"geneName\"";
    int total_reads = 0;
    for(vector<string>::iterator it=MapResultFiles.begin();it!=MapResultFiles.end();it++){
       out<<"\t"<<"\""<<get_file_name(*it)<<"(reads uniquely mapped to gene)"<<"\"";
       out<<"\t"<<"\""<<get_file_name(*it)<<"(reads mapped to gene)"<<"\"";
       out<<"\t"<<"\""<<get_file_name(*it)<<"(all reads)"<<"\"";
       total_reads += readsCount[*it];
    }
    //out<<"total reads in all files"<<endl;
    out<<endl;
    for(vector<string>::iterator it=geneNames.begin();it!=geneNames.end();it++){
       out<<*it;
       for(vector<string>::iterator it2=MapResultFiles.begin();it2!=MapResultFiles.end();it2++){
           out<<"\t"<<genesExp[*it2][*it].first<<"\t"<<genesExp[*it2][*it].second;
           out<<"\t"<<readsCount[*it2];
       }
       //out<<total_reads<<endl;
       out<<endl;
    }
}
void printResult3(string outputFile, vector<string> MapResultFiles, map<string, map<string, pair<int,int> > > & genesExp,
                  vector<string> geneNames, map<string, int> & readsCount, map<string, int> geneAveLen){
    ofstream out(outputFile.c_str());
    if(!out) {
      Rprintf("cannot open output file %s \n", outputFile.c_str());
      return;
    }
    out<<"\"geneName\"";
    int total_reads = 0;
    for(vector<string>::iterator it=MapResultFiles.begin();it!=MapResultFiles.end();it++){
       out<<"\t"<<"\""<<get_file_name(*it)<<"(raw counts)"<<"\"";
       out<<"\t"<<"\""<<get_file_name(*it)<<"(RPKM)"<<"\"";
       out<<"\t"<<"\""<<get_file_name(*it)<<"(all reads)"<<"\"";
       total_reads += readsCount[*it];
    }
    out<<"\t"<<"\"gene length (average of all possible isoform's length)\""<<endl;
    for(vector<string>::iterator it=geneNames.begin();it!=geneNames.end();it++){
       out<<*it;
       for(vector<string>::iterator it2=MapResultFiles.begin();it2!=MapResultFiles.end();it2++){
           out<<"\t"<<genesExp[*it2][*it].first;
           out<<"\t"<<(((((double)genesExp[*it2][*it].first)*1000)/readsCount[*it2])*1000000)/geneAveLen[*it];
           out<<"\t"<<readsCount[*it2];
       }
       out<<"\t"<<geneAveLen[*it];
       out<<endl;
    }
}
int getCounts(string refFlatFile, vector<string> MapResultFiles, string outputFile, string Format, bool needSameStrand){
    vector<Exon> exons;
    map<string,int>  maxExonLen; // chromoson --> int
    vector<string> geneNames;
    map<string, int> readsCount;
    map<string, int> geneAveLen;
    if(getExons(refFlatFile, exons, geneNames, maxExonLen, geneAveLen) < 0){
    	 Rprintf("There is something wrong!\n");
       Rprintf("Please check %s!\n", refFlatFile.c_str());
       return -1;
    }
    map<string, map<string, pair<int,int> > > genesExp;
    for(vector<string>::iterator it=MapResultFiles.begin();it!=MapResultFiles.end();it++){
       map<string, pair<int,int> > tmp;
       genesExp.insert(map<string, map<string, pair<int,int> > >::value_type(*it,tmp));
       if(getCountsForOneSample(exons, maxExonLen, geneNames, *it, Format, needSameStrand, genesExp[*it], readsCount) < 0){
          overlap_percent = 1;
          addtionCol = false;
          return -1;
       }
    }
    printResult3(outputFile, MapResultFiles, genesExp, geneNames, readsCount, geneAveLen);
    overlap_percent = 1;
    addtionCol = false;
    return 0;
}
int getCountsForOneSample(vector<Exon> exons, map<string,int> maxExonLen, vector<string> geneNames,
                          string MapResultFile, string Format, bool needSameStrand,
                          map<string,pair<int,int> > & genesExp, map<string, int> & readsCount){
    //vector<Exon> exons;
    //map<string,int>  maxExonLen; // chromoson --> int
    //vector<string> geneNames;
    //getExons(refFlatFile, exons, geneNames, maxExonLen);
    //map<string,pair<int,int> > genesExp;
    string basename = get_file_name(MapResultFile);
    genesExp.clear();
    for(vector<string>::iterator it=geneNames.begin();it!=geneNames.end();it++){
       pair<int,int> tmp(0,0);
       genesExp.insert(map<string,pair<int,int> >::value_type(*it,tmp));
    }

    ifstream in(MapResultFile.c_str());
    if(!in) {
      Rprintf("cannot open input file %s\n", MapResultFile.c_str());
      R_FlushConsole();
      return -1;
    }
    int length=0;
    if(Format.find("eland")!= string::npos){
      length=atoi((Format.substr(5)).c_str());
    }else{
      length=0;
    }
    while(!in.eof()){
      char buffer[10000 + 1];
      in.getline(buffer, 10000);
      if (buffer[strlen(buffer) - 1] == '\r'){
          buffer[strlen(buffer) - 1] = '\0';
      }
      vector<string> blocks;
      string tmp = buffer;
      if(tmp.size()<10){
         continue;
      }
      string2tokens(string(buffer),blocks);
      if(blocks.size()<5){
         continue;
      }
      if((blocks.size()<6)&&(Format=="bed")){
         continue;
      }
      if((blocks.size()<6)&&(addtionCol==true)){
         continue;
      }
      if((blocks.size()<7)&&(Format=="bed")&&(addtionCol==true)){
         continue;
      }
      if(rows%10000==0){
         Rprintf("\rprocessed %d reads (%s)",rows, basename.c_str());
	       R_FlushConsole();
      }
      rows++;
      string chr;
      int begin;
      int end;
      int strand=0;
      int count=1;
      if(Format=="bed"){
            chr=blocks[0];
            begin=atoi(blocks[1].c_str());
            end=atoi(blocks[2].c_str());
            if(blocks[5]=="+"){
               strand=0;
            }else if(blocks[5]=="-"){
               strand=1;
            }else{
               Rprintf("Wrong Format!\n");
	       R_FlushConsole();
               rows=0;
               overlap_percent = 1;
               addtionCol = false;
               return -1;
            }
            if(addtionCol==true){
               count=atoi(blocks[5].c_str());
            }
      }else if(Format.find("eland")!= string::npos){
            chr=blocks[1];
            string chromoson(chr.substr(0,chr.find_last_of(".fa")-3+1));
            chr=chromoson;
            begin=atoi(blocks[2].c_str())-1;
            end=begin+length;
            if(blocks[4]=="F"){
               strand=0;
            }else if(blocks[4]=="R"){
               strand=1;
            }else{
               Rprintf("Wrong Format!\n");
	             R_FlushConsole();
               rows=0;
               overlap_percent = 1;
               addtionCol = false;
               return -1;
            }
            if(addtionCol==true){
               count=atoi(blocks[4].c_str());
            }
       }else{
            Rprintf("Wrong Format!\n");
	          R_FlushConsole();
            rows=0;
            overlap_percent = 1;
            addtionCol = false;
            return -1;
       }
       Read one_read(chr,begin,end,strand,count);
       findOverLapGenes(exons, one_read, genesExp, maxExonLen,needSameStrand);
    }
    Rprintf("\rprocessed %d reads (%s) \n",rows, basename.c_str());
    R_FlushConsole();
    //printResult(outputFile,genesExp,geneNames);
    readsCount.insert(map<string, int>::value_type(MapResultFile, rows));
    rows=0;
    // overlap_percent = 1;
    // addtionCol = false;
    return 0;
}
int getExonAnnotationFile2(char ** refFlatFile_str, char ** exonAnnotationFile)
{
  string refFlatFile = refFlatFile_str[0];
  string outputFile = exonAnnotationFile[0];
  clock_t start = clock();
  Rprintf("Generate annotation file for exons.\n");
  Rprintf("This will take several minutes.\nPlease wait ...\n");
  R_FlushConsole();
  vector<Exon> exons;
  map<string,int>  maxExonLen; // chromoson --> int
  vector<string> geneNames;
  map<string, int> readsCount;
  map<string, int> geneAveLen;
  if(getExons(refFlatFile, exons, geneNames, maxExonLen, geneAveLen) < 0){
     Rprintf("There is something wrong!\n");
     Rprintf("Please check %s!\n", refFlatFile.c_str());
     return -1;
  }
  ofstream out(outputFile.c_str());
  if(!out) {
     Rprintf("cannot open output file %s \n", outputFile.c_str());
     return -1;
  }
  for(vector<Exon>::iterator it=exons.begin();it!=exons.end();it++){
      string strand;
      if(it->strand == 0){
         strand = "+";
      }else{
         strand = "-";
      }
      out<<it->geneName<<"_"<<it->chr<<"_"<<it->start<<"_"<<it->end<<"_"<<it->strand<<"\t";
      out<<it->geneName<<"_"<<it->chr<<"_"<<it->start<<"_"<<it->end<<"_"<<it->strand<<"\t";
      out<<it->chr<<"\t"<<strand<<"\t"<<it->start<<"\t"<<it->end<<"\t"<<it->start<<"\t"<<it->end<<"\t";
      out<<"1\t"<<it->start<<",\t"<<it->end<<",\n";

  }
  clock_t end = clock();
  Rprintf("total %d unique exons\n", exons.size());
  Rprintf("total used %f seconds.\n",((double)end - start)/CLOCKS_PER_SEC);
  R_FlushConsole();
  return 0;
}
extern "C" {
  static R_NativePrimitiveArgType getGeneExp_t[] = {STRSXP, STRSXP, INTSXP, STRSXP, STRSXP, INTSXP, INTSXP, REALSXP};
  static R_NativePrimitiveArgType getExonAnnotationFile_t[] = {STRSXP, STRSXP};
  R_CMethodDef cMethods[] = {
    {".getGeneExp", (DL_FUNC) &getGeneExp, 8, getGeneExp_t},
    {".getExonAnnotationFile", (DL_FUNC) &getExonAnnotationFile, 2, getExonAnnotationFile_t},
    {NULL, NULL, 0}
  };
  void R_init_DEGseq(DllInfo *info){
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  }
}
