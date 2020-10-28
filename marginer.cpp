// MARGINER standalone 15th Oct 2020
#include <iostream>
#include <string>
#include <vector>
#include "algorithms.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "margins.hpp"
#include <fstream>

//uint STATICMGN = 255 * 255 * 255 * 255;

typedef unsigned int UINT;

int SIZE_DBL = sizeof(double), SIZE_INT = sizeof(UINT), SIZE_CHAR=sizeof(char), NUMC = 256;

std::string FILENAME = "cspace.mgn";
std::vector< iCom<UINT> > coms;

void setfile(){
    // set up file and load into coms;
    std::fstream fr;
    fr.open(FILENAME, std::ios::binary | std::ios::in); // read first
    if (fr.fail()){
        std::cout << "Read failed...\nCreating file...";
        fr.close();
        std::fstream fw;
        fw.open(FILENAME, std::ios::binary | std::ios::out);
        UINT z_int = 0;
        double z_flt = 0.0;
        char z_chr = '0';
        for (int i=0; i < NUMC; i++){
            iCom<UINT> cmm = iCom<UINT>("cspace", {1});
            coms.push_back(cmm);
            for (int j=0; j < NUMC; j++){
                fw.write((char * ) &z_chr, SIZE_CHAR);
                fw.write((char * ) &z_int, SIZE_INT);
                fw.write((char * ) &z_int, SIZE_INT);
                fw.write((char * ) &z_flt, SIZE_DBL);
            }
        }
        fw.close();
    }else{
        std::cout << "Found 'cspace.mgn'...\nLoading...";
        char isset;
        UINT loadedints[2]; // mgn . interrupt . ipm
        double loadedfloat;
        for (int i=0; i < NUMC; i++){
            iCom<UINT> cmm = iCom<UINT>("cspace", {1});
            for (int j=0; j < NUMC; j++){
                fr.read((char *) &isset, SIZE_CHAR);
                fr.read((char *) &loadedints, 2 * SIZE_INT);
                fr.read((char *) &loadedfloat, SIZE_DBL);
                if (isset != '0'){
                    cmm.setnames((UINT) j);
                    cmm.setmgns(loadedints[0]);
                    cmm.setinterrupts(loadedints[1]);
                    cmm.setipms(loadedfloat);
                }
            }
            if (i==-2){
                std::cout << "\nLoaded :" << i << "|\n"; atl::printvector(cmm.getnames());
                std::cout << "\n"; atl::printvector(cmm.getmgns());
                std::cout << "\n"; atl::printvector(cmm.getinterrupts());
                std::cout << "\n"; atl::printvector(cmm.getipms());
            }
            coms.push_back(cmm);
        }
    }
}

void writefile(){
    std::fstream fw;
    fw.open(FILENAME, std::ios::binary | std::ios::out); // equivalent to new write
    std::vector<UINT> names, mgns, interrupts;
    std::vector<double> ipms;
    for (int i=0; i < NUMC; i++){
        names = coms[i].getnames();
        mgns = coms[i].getmgns();
        interrupts = coms[i].getinterrupts();
        ipms = coms[i].getipms();
        if (i==-2){
            std::cout << "\nWriting...\n";
            atl::printvector(coms[i].getnames());
            std::cout << "\n"; atl::printvector(coms[i].getmgns());
            std::cout << "\n"; atl::printvector(coms[i].getinterrupts());
            std::cout << "\n"; atl::printvector(coms[i].getipms());
        }
        int k=0;
        UINT name, m_, i_;
        float f_; char isset;
        for (int j=0; j < NUMC; j++){
            if (k == names.size()){
                isset = '0';
                f_ = 0.0;
                m_ = 0;
                i_ = 0;
                goto wrsect;
            }
            name = names[k];
            if (j < name){
                isset = '0';
                f_ = 0.0;
                m_ = 0;
                i_ = 0;
            }else{
                isset = '1';
                f_ = ipms[k];
                m_ = mgns[k];
                i_ = interrupts[k];
                k++;
            }
            wrsect:
            fw.write((char *) &isset, SIZE_CHAR);
            fw.write((char *) &m_, SIZE_INT);
            fw.write((char *) &i_, SIZE_INT);
            fw.write((char *) &f_, SIZE_DBL);
        }
    }
    fw.close();
}

// the std dev :: iCom test
int main(){
    std::string imagename_, imagename;
    std::cout << "--MARGINER (Mulindwa Myne : 15/10/2020)--";
    std::cout << "\nEnter Image Path : ";
    std::getline(std::cin, imagename_);
    for (int i=0; i < imagename_.size(); i++){
        if (imagename_[i] != '"'){
            imagename.push_back(imagename_[i]);
        }
    }
    setfile();
    cv::Mat tst, cch, imgg;
    imgg = cv::imread(imagename, 0);
    cv::Mat blank = cv::Mat(imgg.rows, imgg.cols, CV_8U, cv::Scalar(0));
    std::cout << "\nSTARTED...";
    for (int y=0; y < imgg.rows; y++){
        for (int x=0; x < imgg.cols - 1; x++){
            uchar clr = imgg.at<uchar>(y, x);
            uchar nclr = imgg.at<uchar>(y, x+1);
            coms[(UINT)clr].add(nclr);
        }
    }
    for (int i=0; i<256; i++){
        coms[i].marginalize();
        //std::cout << "\n" << i << " -> "; atl::printvector(coms[i].getfacts());
    }
    std::cout << "\nUPDATING...";
    writefile();
    std::cout << "\ndrawing";
    for (int y=0; y < imgg.rows - 1; y++){
        for (int x=0; x < imgg.cols - 1; x++){
            uchar clr = imgg.at<uchar>(y, x);
            uchar nclr = imgg.at<uchar>(y, x+1);
            uchar yclr = imgg.at<uchar>(y+1, x);
            if (coms[clr].testif(nclr) == false){
                blank.at<uchar>(y, x) = 255;
            }
            if (coms[clr].testif(yclr) == false){
                blank.at<uchar>(y, x) = 255;
            }
        }
    }
    cv::namedWindow("--MARGINER(Myne)--", cv::WINDOW_NORMAL);
    cv::imshow("--MARGINER(Myne)--", blank);
    cv::waitKey();
}
